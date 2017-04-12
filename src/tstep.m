classdef tstep < handle
% This class defines the functions required to advance the geometry
% forward in time.

properties

tstep_order             % time stepping order
dt                      % time step size
num_particles           % number of particles
num_walls               % number of walls
points_per_particle     % points per particle
points_per_wall         % points per wall
Dp                      % Stokes double-layer potential for fiber-fiber interaction
Dw                      % Stokes double-layer potential for wall-wall interaction
Dup                     % Upsampled Stokes double-layer potential matrix for fibers
Duw                     % Upsampled Stokes double-layer potential matrix for walls
N0w                     % N0 matrix to remove rank 1 nullspace
fmm                     % flag for using the FMM
near_singular           % flag for using near-singular integration
gmres_tol               % GMRES tolerance
gmres_max_it            % maximum GMRES iterations
near_structff           % near-singular integration structure (fibre-fibre)
near_structfw           % near-singular integration structure (fibre-wall)
near_structwf           % near-singular integration structure (wall-fibre)
far_field               % background flow
use_precond             % use a block-diagonal preconditioner
precop                  % block-diagonal preconditioner for fibres
precow                  % block-diagonal preconditioner for walls
potp                    % class for fiber layer potentials
potw                    % class for wall layer potentials
profile                 % flag to time certain parts of code
om                      % monitor class
confined                % flag to indicate whether flow is bounded
tau0                    % initial angles of fibres
resolve_collisions      % flag to indicate whether to resolve collisions
prams
minimum_separation
display_solution

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(options, prams, om, geom, walls, tau0)
% o.tstep(options,prams): constructor.  Initialize class.  Take all
% elements of options and prams needed by the time stepper

o.tstep_order = options.tstep_order;        
o.fmm = options.fmm;                  
o.near_singular = options.near_singular;                
o.use_precond = options.use_precond;           
o.profile = options.profile;             
o.confined = options.confined;
o.prams = prams;
o.minimum_separation = prams.minimum_separation;
o.display_solution = options.display_solution;

o.resolve_collisions = options.resolve_collisions;
o.dt = prams.T/prams.number_steps;                  
o.gmres_tol = options.gmres_tol;             
o.gmres_max_it = 100;

o.far_field = @(X) o.bgFlow(X,options); 
o.om = om;
o.tau0 = tau0;

o.num_particles = prams.np;
o.num_walls = prams.nw;
o.points_per_particle = prams.Np;
o.points_per_wall = prams.Nw;

np = o.num_particles;
Np = o.points_per_particle;
Nw = o.points_per_wall;

% CREATE CLASSES TO EVALUATE POTENTIALS ON FIBRES AND WALLS
o.potp = poten(np, om);

if options.confined
    o.potw = poten(Nw, om);
    o.Dw = o.potw.stokesDLmatrix(walls);
else
    o.potw = [];
    o.Dw = [];
end

% CREATE MATRICES FOR FIBRE-FIBRE SELF INTERATIONS AND 
% WALL-WALL SELF INTERACTIONS
if ~isempty(geom.X)
    o.Dp = o.potp.stokesDLmatrix(geom);
else
    o.Dp = [];
end
% CREATE UPSAMPLED MATRICES
% FIBRE-FIBRE
if ~isempty(geom.X)
    Xsou = geom.X; 
    Nup = Np*ceil(sqrt(Np));

    Xup = [interpft(Xsou(1:Np,:),Nup);...
       interpft(Xsou(Np+1:2*Np,:),Nup)];

    geomUp = capsules([],Xup);
    o.Dup = o.potp.stokesDLmatrix(geomUp);
else
    o.Duw = [];
end

% WALL-WALL
if options.confined
    Xsou = walls.X; 
    Nup = Nw*ceil(sqrt(Nw));

    Xup = [interpft(Xsou(1:Nw,:),Nup);...
           interpft(Xsou(Nw+1:2*Nw,:),Nup)];

    wallsUp = capsules([],Xup);
    o.Duw = o.potw.stokesDLmatrix(wallsUp);

    % CREATE N0 MATRIX
    o.N0w = o.potw.stokesN0matrix(walls);
end

% CREATE BLOCK-DIAGONAL PRECONDITIONER
if o.use_precond
    
    if o.profile
        tic;
    end
    
    % FIBEE-FIBRE PRECONDITIONER
    o.precop.L = zeros(2*Np+3,2*Np+3,np);
    o.precop.U = zeros(2*Np+3,2*Np+3,np);
    for k = 1:np
        [o.precop.L(:,:,k),o.precop.U(:,:,k)] =...
            lu([-1/2*eye(2*Np)+o.Dp(:,:,k) ...
            [-ones(Np,1);zeros(Np,1)] ...
            [zeros(Np,1);-ones(Np,1)] ...
            [geom.X(end/2+1:end,k)-geom.center(2,k); ...
                -geom.X(1:end/2,k)+geom.center(1,k)];...
            [-geom.sa(:,k)'*2*pi/Np zeros(1,Np) 0 0 0];
            [zeros(1,Np) -geom.sa(:,k)'*2*pi/np 0 0 0];
            [(-geom.X(end/2+1:end,k)'+geom.center(2,k)).*geom.sa(:,k)'*2*pi/Np ...
            (geom.X(1:end/2,k)'-geom.center(1,k)).*geom.sa(:,k)'*2*pi/Np 0 0 0]]);
    end
    
    if o.confined
        
        % WALL-WALL PRECONDITIONER
        o.precow.L = zeros(2*Nw + 3,2*Nw + 3,Nw);
        o.precow.U = zeros(2*Nw + 3,2*Nw + 3,Nw);
        
        oc = curve;
        sa = walls.sa;
        [x,y] = oc.getXY(walls.X);
        [cx,cy] = oc.getXY(walls.center);
        
        for k = 1:Nw
            
            if k == 1 % first wall does not need Rotlets and Stokeslets
                [o.precow.L(1:2*Nw,1:2*Nw,k),o.precow.U(1:2*Nw,1:2*Nw,k)] =...
                    lu(-1/2*eye(2*Nw)+o.Dw(:,:,k)+o.N0w(:,:,k));
                
            else
                r = [x(:,k) - cx(k), y(:,k) - cy(k)];
                rho2 = (x(:,k) - cx(k)).^2 + (y(:,k) - cy(k)).^2;
                
                col_stokes1 = [-0.5*log(rho2) + r(:,1).*r(:,1)./rho2; ...
                                r(:,2).*r(:,1)./rho2]/(4*pi);
                col_stokes2 = [r(:,2).*r(:,1)./rho2; -0.5*log(rho2) +...
                                r(:,2).*r(:,2)./rho2]/(4*pi);
                col_rot = [r(:,2)./rho2; -r(:,1)./rho2];
                
                int_stokes = 2*pi*sa(:,k)'/Nw;
                int_rot = 2*pi*[(y(:,k).*sa(:,k))', -(x(:,k).*sa(:,k))']/Nw;
                
                [o.precow.L(:,:,k),o.precow.U(:,:,k)] =...
                    lu([-1/2*eye(2*Nw)+o.Dw(:,:,k), col_stokes1, ...
                    col_stokes2, col_rot;...
                    int_stokes, zeros(1,Nw), -2*pi, 0, 0;...
                    zeros(1,Nw), int_stokes, 0, -2*pi, 0;...
                    int_rot, 0, 0, -2*pi]);
            end
        end
    end
    
    if o.profile
        o.om.writeMessage(['Building preconditioner ... ', num2str(toc)]);
    end    
end

end % constructor: tstep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xc_new,tau_new,etaP,etaW,uP,omega,forceW,torqueW,forceP,....
                torqueP,iter,iflag,res] = timeStep(o,geom,walls,xc,tau,...
                uP_m1,omega_m1,first_step)
% timeStep(geom, walls, xc, tau) takes a single time step and returns the
% angles and velocities of the particles at the next time step as well as
% the density functions on the particles and walls, the translational and
% rotational velocities of the particles, and the net force and torque on
% all particles and walls

np = o.num_particles;
Np = o.points_per_particle;

% ASSEMBLE RHS WITH NO FORCE AND TORQUE ON PARTICLES
forceP = zeros(2,np);
torqueP = zeros(1,np);

rhs = o.assembleRHS(geom, walls, forceP, torqueP, false);

% CREATE NEAR SINGULAR INTEGRATION STRUCTURES
if o.near_singular
    if o.profile
        tic;
    end
    
    o.near_structff = geom.getZone([],1);
    
    if o.confined
        [~,o.near_structfw] = geom.getZone(walls,2);
        [~,o.near_structwf] = walls.getZone(geom,2);
    end
    
    if o.profile
        o.om.writeMessage(['getZone ... ', num2str(toc)]);
    end
end

% ROTATE FIBRE-FIBRE DLP AND FIBRE-FIBRE PRECONDITIONER
Nup = Np*ceil(sqrt(Np));
dtau = tau - o.tau0;
o.tau0 = tau;

% create diagonal blocks without rotation
% Xup = [interpft(geom.X(1:N,:),Nup);...
%        interpft(geom.X(N+1:2*N,:),Nup)];
% 
% geomUp = capsules([],Xup);
% o.Dupf = o.opFibers.stokesDLmatrix(geomUp);
% 
% o.Df = o.opFibers.stokesDLmatrix(geom);
% o.Dupf = o.opFibers.stokesDLmatrix(geomUp);

for i = 1:np

    R = spdiags([sin(dtau(i))*ones(2*Np,1), cos(dtau(i))*ones(2*Np,1)...
                -sin(dtau(i))*ones(2*Np,1)], [-Np, 0, Np], zeros(2*Np, 2*Np));
    
    Rup = spdiags([sin(dtau(i))*ones(2*Nup,1), cos(dtau(i))*ones(2*Nup,1)...
                -sin(dtau(i))*ones(2*Nup,1)], [-Nup, 0, Nup], zeros(2*Nup, 2*Nup));  
            
    o.Dp(:,:,i) = R*o.Dp(:,:,i)*R';
    o.Dup(:,:,i) = Rup*o.Dup(:,:,i)*Rup';
    
    if o.use_precond
        
        [o.precop.L(:,:,i),o.precop.U(:,:,i)] =...
            lu([-1/2*eye(2*Np)+o.Dp(:,:,i) ...
            [-ones(Np,1);zeros(Np,1)] ...
            [zeros(Np,1);-ones(Np,1)] ...
            [geom.X(end/2+1:end,i)-geom.center(2,i);...
                    -geom.X(1:end/2,i)+geom.center(1,i)];...
            [-geom.sa(:,i)'*2*pi/Np zeros(1,Np) 0 0 0];
            [zeros(1,Np) -geom.sa(:,i)'*2*pi/Np 0 0 0];
            [(-geom.X(end/2+1:end,i)'+geom.center(2,i)).*geom.sa(:,i)'*2*pi/Np ...
            (geom.X(1:end/2,i)'-geom.center(1,i)).*geom.sa(:,i)'*2*pi/Np 0 0 0]]);
    end    
end

% SOLVE SYSTEM USING GMRES TO GET CANDIDATE TIME STEP
if o.use_precond
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom,walls,true),rhs,[],...
      o.gmres_tol,o.gmres_max_it,@o.preconditionerBD,[], rhs);
else
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom,walls,true),rhs,[],...
      o.gmres_tol, o.gmres_max_it);
end

iter = I(2);

% UNPACK SOLUTION
[etaP, etaW, uP, omega, forceW, torqueW] = o.unpackSolution(Xn);

% CREATE CANDIDATE CONFIGURATION 
if first_step % FORWARD EULER
    xc_new = xc + o.dt*uP;
    tau_new = tau + o.dt*omega;
else %ADAMS_BASHFORTH
    xc_new = xc + (3/2)*o.dt*uP - (1/2)*o.dt*uP_m1;
    tau_new = tau + (3/2)*o.dt*omega - (1/2)*o.dt*omega_m1;
end

geomProv = capsules(o.prams, xc_new, tau_new);
    
% RESOLVE COLLISIONS
if o.resolve_collisions
    
    % REASSEMBLE OLD CONFIGURATION
    geomOld = capsules(o.prams,xc, tau);
    
    solution.xc_old = xc;
    solution.tau_old = tau;
    solution.xc_new = xc_new;
    solution.tau_new = tau_new;
    solution.etaP = etaP;
    solution.etaW = etaW;
    solution.uP = uP;
    solution.omega = omega;
    solution.forceW = forceW;
    solution.torqueW = torqueW;
    solution.forceP = 0;
    solution.torqueP = 0;
    
    solution =  o.resolveCollisions(geomOld, geomProv, walls, solution,...
                    uP_m1, omega_m1, first_step);
    
    xc_new = solution.xc_new;
    tau_new = solution.tau_new;
    etaP = solution.etaP;
    etaW = solution.etaW;
    uP = solution.uP;
    omega = solution.omega;
    forceW = solution.forceW;
    torqueW = solution.torqueW;
    forceP = solution.forceP;
    torqueP = solution.torqueP;
end

if o.display_solution
    geomNew = capsules(o.prams, xc_new, tau_new);
    fill(geomNew.X(1:end/2,:),geomNew.X(end/2+1:end,:),'k');
    axis equal
    drawnow
end

end % timeStep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tx = timeMatVec(o,Xn,geom,walls,include_far_field)
% Tx = timeMatVec(Xn,geom) does a matvec for GMRES 
if o.profile
    tMatvec = tic;
end

np = o.num_particles;
nw = o.num_walls;
Np = o.points_per_particle;
Nw = o.points_per_wall;

pot_particles = o.potp;
pot_walls = o.potw;

% PREALLOCATE OUTPUT VECTORS
velParticles = zeros(2*Np,np);
velWalls = zeros(2*Nw,nw);
forceP = zeros(2,np);
torqueP = zeros(np,1);

% UNPACK Xn
[etaP, etaW, uP, omega, forceW, torqueW] = unpackSolution(o, Xn);

% CALCULATE VELOCITIES ON PARTICLES AND WALLS

% ADD JUMP IN DLP

velParticles = velParticles - 1/2*etaP;

if o.confined
   velWalls = velWalls - 1/2*etaW;
end 

% ADD SELF CONTRIBUTION
velParticles = velParticles + ...
    pot_particles.exactStokesDLdiag(geom, o.Dp, etaP);

if o.confined
    velWalls = velWalls + pot_walls.exactStokesDLdiag(walls, o.Dw, etaW);
end

% START OF SOURCE == PARTICLES
% START OF TARGET == PARTICLES
if o.fmm
  kernel = @pot_particles.exactStokesDLfmm;
else
  kernel = @pot_particles.exactStokesDL;
end

kernelDirect = @pot_particles.exactStokesDL;

if o.near_singular 
  DLP = @(X) pot_particles.exactStokesDLdiag(geom,o.Dp,X) - 1/2*X;
  pp_dlp = pot_particles.nearSingInt(geom, etaP, DLP, o.Dup, ...
      o.near_structff ,kernel, kernelDirect, geom, true, false);
else
  pp_dlp = kernel(geom, etaP, o.Dp);
end
% END OF TARGET == PARTICLES

% START OF TARGET == WALLS 
if o.confined && np > 0
   
   if o.fmm
       kernel = @pot_walls.exactStokesDLfmm;       
   else
       kernel = @pot_walls.exactStokesDL;
   end
   
   kernelDirect = @pot_walls.exactStokesDL;
   
   if o.near_singular
       DLP = @(X) pot_walls.exactStokesDLdiag(geom, o.Dp, X) - 1/2*X;
       wp_dlp = pot_walls.nearSingInt(geom, etaP, DLP, o.Dup, ...
           o.near_structfw, kernel, kernelDirect, walls, false, false);
   else
       wp_dlp = kernel(geom, etaP);
   end
else
    wp_dlp = zeros(2*Nw,nw);
end
% END OF TARGET == WALLS

% START OF SOURCE == WALLS
% START OF TARGET == PARTICLES
if o.confined && np > 0
   
   if o.fmm
        kernel = @pot_particles.exactStokesDLfmm;       
   else
       kernel = @pot_particles.exactStokesDL;
   end
   
   kernelDirect = @pot_particles.exactStokesDL;
   
   if o.near_singular
       DLP = @(X) pot_particles.exactStokesDLdiag(walls, o.Dw, X) - 1/2*X;
       pw_dlp = pot_particles.nearSingInt(walls, etaW, DLP, o.Duw, ...
           o.near_structwf, kernel, kernelDirect, geom, false, false);
   else
       pw_dlp = kernel(walls, etaW);
   end
else
    pw_dlp = zeros(2*Np,np);
end
% END OF TARGET == PARTICLES

% START OF TARGET == WALLS
if o.confined
    
    if o.fmm
        kernel = @pot_walls.exactStokesDLfmm;
    else
        kernel = @pot_walls.exactStokesDL;
    end
    
    ww_dlp = kernel(walls, etaW, o.Dw);
    
else
    ww_dlp = zeros(2*Nw,nw);
end
% END OF TARGET == WALLS
% END OF SOURCE == WALLS

if (o.confined && nw > 1)
    % START SOURCE == ROTLETS AND STOKESLETS
    % START TARGETS == PARTICLES
    p_sr = 0;
    for k = 2:nw % loop over all walls, except outer wall
        p_sr = p_sr + ...
            o.RSlets(geom.X, walls.center(:,k), forceW(:,k-1), torqueW(k-1));
    end 
    
    % END TARGETS == PARTICLES
    
    % START TARGETS == WALLS
    w_sr = zeros(2*Nw,nw);
    for k = 2:nw % loop over all walls, except outer wall
        w_sr = w_sr +...
            o.RSlets(walls.X, walls.center(:,k), forceW(:,k-1), torqueW(k-1));
    end   
    
    % END TARGETS == WALLS
    % END SOURCE == ROTLETS
    
    % EVALUATE ROTLET AND STOKESLET EQUATIONS
    valLets = o.letsIntegrals([forceW;torqueW], etaW, walls);
else
    p_sr = 0;
    w_sr = 0;
    valLets = [];
end

% EVALUATE TOTAL VELOCITY ON PARTICLES

% ADD CONTRIBUTIONS FROM OTHER BODIES
if include_far_field
    velParticles = velParticles + p_sr + pp_dlp + pw_dlp;
end

% SUBTRACT VELOCITY ON SURFACE
for k = 1:np
  velParticles(1:end/2,k) = velParticles(1:end/2,k) - uP(1,k);
  velParticles(end/2+1:end,k) = velParticles(end/2+1:end,k) - uP(2,k);
end

for k = 1:np
  velParticles(1:end/2,k) = velParticles(1:end/2,k) ...
                - (geom.X(end/2+1:end,k) - geom.center(2,k))*omega(k);
  velParticles(end/2+1:end,k) = velParticles(end/2+1:end,k)...
                + (geom.X(1:end/2,k) - geom.center(1,k))*omega(k); 
end

% EVALUATE VELOCITY ON WALLS
if o.confined
    velWalls = velWalls + ww_dlp + wp_dlp + w_sr;
    velWalls(:,1) = velWalls(:,1)...
                    + pot_walls.exactStokesN0diag(walls, o.N0w, etaW(:,1));
end

% EVALUTATE FORCES ON PARTICLES
for k = 1:np
  forceP(1,k) = sum(etaP(1:Np,k).*geom.sa(:,k))*2*pi/Np;
  forceP(2,k) = sum(etaP(Np+1:2*Np,k).*geom.sa(:,k))*2*pi/Np;
end

% EVALUATE TORQUES ON PARTICLES
for k = 1:np
  torqueP(k) = sum(((geom.X(Np+1:2*Np,k)-geom.center(2,k)).*etaP(1:Np,k) - ...
                     ((geom.X(1:Np,k)-geom.center(1,k)).*etaP(Np+1:2*Np,k))).*...
                     geom.sa(:,k))*2*pi/Np;
end

% CONSTRUCT OUTPUT VECTOR
Tx = [velParticles(:); velWalls(:); forceP(:); torqueP(:); valLets(:)];

if o.profile
    o.om.writeMessage(['Matvec assembly completed in ', ....
                            num2str(toc(tMatvec)), ' seconds']);
end

end % timeMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhs = assembleRHS(o, geom, walls, forceP, torqueP, preprocess)
% Assemble right hand side

np = o.num_particles;
nw = o.num_walls;
Np = o.points_per_particle;
Nw = o.points_per_wall;

if o.confined
    ff = o.far_field(walls.X);
    rhs = [zeros(2*Np*np,1); ff(:); forceP(:); torqueP' ;zeros(3*(Nw-1),1)];
else
    
    ff = o.far_field(geom.X);
    rhs = [-ff(:); forceP(:); torqueP'];
end

if o.resolve_collisions
    
    % ADD IN CONTRIBUTIONS TO VELOCITIES FROM PARTICLE ROTLETS AND STOKESLETS    
    for k = 1:np
        % CONTRIBUTION TO PARTICLE VELOCITY
        v = o.RSlets(geom.X, geom.center(:,k),forceP(:,k),torqueP(k));
        
        if preprocess
            start = 2*Np*(k-1)+1;
            rhs(start:start+2*Np-1) = rhs(start:start+2*Np-1) - v(:,k);
        else
            rhs(1:2*Np*np)= rhs(1:2*Np*np) - v(:);
        end
        
        if o.confined
            % CONTRIBUTION TO WALL VELOCITY
            v = o.RSlets(walls.X, geom.center(:,k),forceP(:,k),torqueP(k));
            rhs(2*Np*np+1:2*Np*np+2*Nw*nw) = ...
                                rhs(2*Np*np+1:2*Np*np+2*Nw*nw) - v(:);
        end
    end
end

end % assembleRHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [etaP, etaW, u, omega, lambda, xi] = unpackSolution(o, Xn)
% Unpack solution vector that comes from GMRES call

np = o.num_particles;
nw = o.num_walls;
Np = o.points_per_particle;
Nw = o.points_per_wall;

% REORGANIZE COLUMN VECTOR INTO MATRIX
% EXTRACT DENSITY FUNCITONS ON PARTICLES AND WALLS
% each column of etaP corresponds to the density function of a rigid body
etaP = zeros(2*Np,np);
for k = 1:np
  etaP(:,k) = Xn((k-1)*2*Np+1:k*2*Np);
end

% each column of etaW corresponds to the density function of a solid wall
if o.confined
    etaW = zeros(2*Nw, nw);
    for k = 1:nw
       etaW(:,k) = Xn(2*Np*np+1+(k-1)*2*Nw:2*Np*np+k*2*Nw);
    end
else
    etaW = 0;
end

% EXTRACT TRANSLATIONAL AND ROTATIONAL VELOCITIES
u = zeros(2,np);
for k = 1:np
  u(:,k) = Xn(2*Np*np+2*Nw*nw+(k-1)*2+1:2*Np*np+2*Nw*nw+k*2);
end

omega = zeros(1,np);
for k = 1:np
  omega(k) = Xn(2*Np*np+2*Nw*nw+2*np+k);
end

if o.confined
    lambda = zeros(2,nw-1);
    for k = 1:nw-1
       lambda(:,k) = ...
           Xn(2*Np*np+2*Nw*nw+3*np+1+2*(k-1): 2*Np*np+2*Nw*nw+3*np+2*k);
    end

    xi = zeros(1,nw-1);
    for k = 1:nw-1
       lambda(k) = Xn(2*Np*np+2*Nw*nw+3*np+2*(nw-1)+k);
    end
else
    lambda = 0;
    xi = 0;
end

end % unpackSolution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pz = preconditionerBD(o,z)
% apply the block-diagonal preconditioner whose LU factorization is
% precomputed and stored

Np = o.points_per_particle;
np = o.num_particles;
Nw = o.points_per_wall;
nw = o.num_walls;

Pz = z;

% APPLY FIBRE-FIBRE PRECONDITIONER
for k = 1:np
  etaStart = (k-1)*2*Np + 1;
  endEnd = etaStart + 2*Np - 1;
  uStart = 2*Np*np + 2*Nw*nw + (k-1)*2 + 1;
  uEnd = uStart + 1;
  omegaStart = 2*Np*np + 2*Nw*nw + np*2 + (k-1) + 1;
  omegaEnd = omegaStart;
  zblock = o.precop.U(:,:,k)\(o.precop.L(:,:,k)\...
      [z(etaStart:endEnd);z(uStart:uEnd);z(omegaStart:omegaEnd)]);

  Pz(etaStart:endEnd) = zblock(1:2*Np);
  Pz(uStart:uEnd) = zblock(2*Np+1:2*Np+2);
  Pz(omegaStart:omegaEnd) = zblock(2*Np+3:2*Np+3);
end

% APPLY WALL-WALL PRECONDITIONER
for k = 1:nw
    
    if k == 1 %No stokeslets or rotlets
       xiStart = 2*Np*np +  1;
       xiEnd = xiStart + 2*Nw - 1;
       
       zblock = o.precow.U(1:2*Nw,1:2*Nw,1)\(o.precow.L(1:2*Nw,1:2*Nw,1)\...
                            z(xiStart:xiEnd));
       Pz(xiStart:xiEnd) = zblock;
       
    else
        xiStart = 2*Np*np + 2*(k-1)*Nw + 1;
        xiEnd = xiStart + 2*Nw - 1;
        
        stokesletStart = 2*Np*np+2*Nw*nw+3*np+2*(k-2)+1;
        stokesletEnd = stokesletStart + 1;
        rotletStart = 2*Np*np+2*Nw*nw+3*np+2*(nw-1)+(k-1);
        rotletEnd = rotletStart;
        
        zblock = o.precow.U(:,:,k)\(o.precow.L(:,:,k)\...
                [z(xiStart:xiEnd);z(stokesletStart:stokesletEnd);...
                z(rotletStart:rotletEnd)]);

        Pz(xiStart:xiEnd) = zblock(1:2*Nw);
        Pz(stokesletStart:stokesletEnd) = zblock(2*Nw+1:2*Nw+2);
        Pz(rotletStart:rotletEnd) = zblock(2*Nw+3:2*Nw+3);
    end
end

end % preconditioner

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = RSlets(~,X,center,stokeslet,rotlet)
% vel = RSlets(o,X,center,stokeslet,rotlet) evaluates the velocity due
% to the stokeslet and rotlet terms at target points X.  
% Center of the rotlet and stokeslet is contained in center

oc = curve;
[x,y] = oc.getXY(X);
% set of points where we are evaluating the velocity
[cx,cy] = oc.getXY(center);
% the center of the rotlet/stokeslet terms

rho2 = (x-cx).^2 + (y-cy).^2;
% distance squared

LogTerm = -0.5*log(rho2)*stokeslet(1);
rorTerm = 1./rho2.*((x-cx).*(x-cx)*stokeslet(1) + ...
    (x-cx).*(y-cy)*stokeslet(2));
RotTerm = (y-cy)./rho2*rotlet;
velx = 1/4/pi*(LogTerm + rorTerm) + RotTerm;
% x component of velocity due to the stokeslet and rotlet

LogTerm = -0.5*log(rho2)*stokeslet(2);
rorTerm = 1./rho2.*((y-cy).*(x-cx)*stokeslet(1) + ...
    (y-cy).*(y-cy)*stokeslet(2));
RotTerm = -(x-cx)./rho2*rotlet;
vely = 1/4/pi*(LogTerm + rorTerm) + RotTerm;
% y component of velocity due to the stokeslet and rotlet

vel = [velx;vely];
% velocity

end % RSlets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = letsIntegrals(~,otlets,eta,walls)
% z = letsIntegrals(stokeslet,rotlet,etaM,walls) integrates the density
% function to enforce constraints on stokeslets and rotlets

Nw = o.points_per_wall;
nw = o.num_walls;

z = zeros(3*(nw-1),1);

for k = 2:nw
  stokeslet = otlets(3*(k-2)+1:3*(k-2)+2);
  % two stokeslet terms per inner boundary
  rotlet = otlets(3*(k-1));
  % one rotlet term per inner boundary
  ind = 3*(k-2)+1;
  z(ind) = -2*pi*stokeslet(1) + ...
    sum(eta(1:Nw,k).*walls.sa(:,k))*2*pi/Nw;
  % integral of density function dotted with [1;0]
  % is one stokeslet
  z(ind+1) = -2*pi*stokeslet(2) + ...
    sum(eta(Nw+1:2*Nw,k).*walls.sa(:,k))*2*pi/Nw;
  % integral of density fuction dotted with [0;1]
  % is the other stokeslet
  z(ind+2) = -2*pi*rotlet + sum(...
    ((walls.X(Nw+1:2*Nw,k)).*eta(1:Nw,k) - ...
    (walls.X(1:Nw,k)).*eta(Nw+1:2*Nw,k)).*...
    walls.sa(:,k))*2*pi/Nw;
  % integral of density function dotted with (y,-x)
  % is the rotlet
end % k

end % letsIntegrals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solution = resolveCollisions(o,geomOld, geomProv, walls, ...
                        solution, uP_m1, omega_m1, first_step)

Np = o.points_per_particle;
np = o.num_particles;
Nw = o.points_per_wall;
nw = o.num_walls;

oc = curve;
cellSize = 0;
if np
  [~,length] = oc.geomProp(geomOld.X);
  edgelength = length/Np;
  cellSize = max(cellSize,max(edgelength));
end

if o.confined
  [~,length] = oc.geomProp(walls.X);
  walllength = max(length/Nw);
  wallupSamp = ceil(walllength/cellSize);
else
  wallupSamp = 1;
end

upSampleFactor = 1;
nexten = 0;
c_tol = 1e-12;
minSep = o.minimum_separation;
maxIter = 1000;

% CHECK FOR COLLISION
[Ns,totalnp,Xstart,Xend,totalPts] = o.preColCheck(geomOld.X,geomProv.X,...
                                        walls,wallupSamp);
[vgrad, iv, ids, vols] = getCollision(Ns, totalnp, Xstart, Xend, minSep, ...
        maxIter, totalPts, c_tol, np, nw, Np*upSampleFactor, ...
        Nw*wallupSamp, nexten, max(cellSize,minSep));
    
o.om.writeMessage(['ivolume: ' num2str(iv)],'%s\n');

fc_tot = zeros(2*Np*np,1);
colCount = 0;

% RESOLVE COLLISION
while(iv<0)

    colCount = colCount + 1;  

    vgrad = vgrad(1:2*Np*np);

    disp(max(abs(vgrad(1:Np))));
    ids = ids(1:2*Np*np);
    vols = vols(1:2*Np*np);

    vgrad(1:2*Np*np) = o.adjustNormal(vgrad(1:2*Np*np),Np,np,geomOld,...
                                edgelength,colCount);
    [A,~,ivs,~,jacoSmooth] = o.preprocessRigid(vgrad,ids,vols,geomOld,walls);

    % CALCULATE FORCE AND TORQUES ON PARTICLES
    [fc, ~] = o.getColForce(A,ivs/o.dt,ivs*0,jacoSmooth);
    fc_tot = fc_tot + fc;

    fc_tot_tmp = reshape(fc_tot, 2*Np, np);
    [forceP,torqueP] = ...
            o.getRS(geomOld.X,geomOld.sa,geomOld.center,fc_tot_tmp);

    % ASSEMBLE NEW RHS WITH CONTACT FORCES
    rhs = o.assembleRHS(geomOld, walls, forceP, torqueP, true);

    % SOLVE SYSTEM
    if o.use_precond
      Xn = gmres(@(X) o.timeMatVec(X,geomOld,walls,true),rhs,[],o.gmres_tol,...
          o.gmres_max_it,@o.preconditionerBD,[]);
    else
      Xn = gmres(@(X) o.timeMatVec(X,geomOld,walls,true),rhs,[],o.gmres_tol,...
          o.gmres_max_it);
    end 

    % UNPACK SOLUTION
    [etaP, etaW, uP, omega, forceW, torqueW] = o.unpackSolution(Xn);

    % UPDATE PARTICLE POSTIONS AND ANGLES 
    if first_step % FORWARD EULER
        xc_new = solution.xc_old + o.dt*uP;
        tau_new = solution.tau_old + o.dt*omega;
    else %ADAMS_BASHFORTH
        xc_new = solution.xc_old + (3/2)*o.dt*uP - (1/2)*o.dt*uP_m1;
        tau_new = solution.tau_old + (3/2)*o.dt*omega - (1/2)*o.dt*omega_m1;
    end

    geomProv = capsules(o.prams, xc_new, tau_new);

    % CHECK AGAIN FOR COLLISIONS
    [Ns,totalnp,Xstart,Xend,totalPts] = o.preColCheck(geomOld.X,geomProv.X,...
      walls,wallupSamp);

    [vgrad, iv, ids, vols] = getCollision(Ns, totalnp, Xstart, Xend, minSep, ...
      maxIter, totalPts, c_tol, np, nw, Np*upSampleFactor, ...
      Nw*wallupSamp, nexten, max(cellSize,minSep));

    solution.xc_new = xc_new;
    solution.tau_new = tau_new;
    solution.etaP = etaP;
    solution.etaW = etaW;
    solution.uP = uP;
    solution.omega = omega;
    solution.forceW = forceW;
    solution.torqueW = torqueW;
    solution.forceP = forceP;
    solution.torqueP = torqueP;

    o.om.writeMessage(['ivolume: ' num2str(iv)],'%s\n');
end


end % resolveCollision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,jaco,ivs,listnv,jacoSmooth] = preprocessRigid(o,vgrad,ids,...
                    vols,geom,walls)
                
Np = o.points_per_particle;
np = o.num_particles;

nivs = max(ids);
A = zeros(nivs,nivs);
ivs = zeros(nivs,1);
vtoiv = zeros(np,nivs);

jI = [];
jJ = [];
jV = [];
listnv = [];

for i = 1:2*Np*np
    if(ids(i)~=0)
        k = ceil(i/(2*Np));  % determine particle number
        listnv = [listnv;k]; % keep running track of all particle involved
        jI = [ids(i);jI];    % keep track of all interference volumes
        jJ = [i;jJ];         % keep track of all nodes 
        jV = [vgrad(i);jV];  % keep track of all volume gradients
        vtoiv(k,ids(i)) = 1; 
        ivs(ids(i)) = vols(i);
    end
end

listnv = unique(listnv);
ivs(ivs>-1e-12) = -1e-12;
% sparse matrix of size nivs by total number of points on particles
jaco = sparse(jI,jJ,jV,nivs,2*Np*np); 
jacoSmooth = jaco*1;

% LOOP OVER INTERSECTION VOLUMES
for i = 1:nivs
    
    % FIND BODIES IN THIS VOLUME
    S = find(vtoiv(:,i)~=0); 
    
    f = jaco(i,1:2*Np*np)';
    f = reshape(f, 2*Np, np);

    [forceP,torqueP] = o.getRS(geom.X,geom.sa,geom.center,full(f));
    
    % ASSEMBLE NEW RHS WITH CONTACT FORCES
    rhs = o.assembleRHS(geom, walls, forceP, torqueP, true);
    ff = o.far_field(geom.X);
    rhs(1:2*Np*np) = rhs(1:2*Np*np) + ff(:);

    % SOLVE SYSTEM
    if o.use_precond
      Xn = gmres(@(X) o.timeMatVec(X,geom,walls,false),rhs,[],o.gmres_tol,...
          o.gmres_max_it,@o.preconditionerBD,[]);
    else
      Xn = gmres(@(X) o.timeMatVec(X,geom,walls,false),rhs,[],o.gmres_tol,...
          o.gmres_max_it);
    end 
    
    % COMPUTE VELOCITY OF EACH POINT ON ALL RIGID PARTICLES
    [~, ~, uP, omega, ~, ~] = unpackSolution(o, Xn);
    
    b = zeros(2*Np,np);
    for k = S'
        b(:,k) = [uP(1,k)*ones(Np,1); uP(2,k)*ones(Np,1)] + ... %translational
            omega(k)*[geom.X(Np+1:end,k) - geom.center(2,k); ...% + rotational
            -(geom.X(1:Np,k) - geom.center(1,k))];
    end

    b = b(:);
    
    SS = find(vtoiv(k,:)~=0);
    for l = 1:numel(SS)
        A(SS(l),i) = A(SS(l),i) + dot(jaco(SS(l),:),b);
    end
end

end % preprocessRigid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vgrad] = adjustNormal(~,vgrad,N,n,geom,edgelength,colCount)
% use the normal of Xn or X?
vgrad = reshape(vgrad,2*N,n);

for k = 1:n
  for i = 1:N
    if(vgrad(i,k)==0 && vgrad(N+i,k)==0)
      continue;
    end
    
    % in rare case if space-time gradient is too small, 
    % collision resolving may get stuck
    % use surface normal instead
    if(colCount > 20)
      bnnorm = edgelength(k);
      
      tangent = geom.xt;
      oc = curve;
      [tanx,tany] = oc.getXY(tangent);
      nx = tany;
      ny = -tanx;
      normal = [nx;ny];
      
      vgrad(i,k) = normal(i,k)*bnnorm;
      vgrad(N+i,k) = normal(i+N,k)*bnnorm;
    end
  end
end

vgrad = vgrad(:);
end % adjustNormal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokeslet, rotlet] = getRS(~,X,sa,cm,eta)
  N = size(X,1)/2;
  n = size(X,2);
  
  stokeslet = zeros(2,n);
  
  stokeslet(1,:) = sum(eta(1:N,:).*sa)*2*pi/N;
  stokeslet(2,:) = sum(eta(N+1:end,:).*sa)*2*pi/N;
  
  rotlet = sum(((X(N+1:end,:)-repmat(cm(2,:),[N,1])).*eta(1:N,:) - ...
    (X(1:N,:)-repmat(cm(1,:),[N,1])).*eta(N+1:end,:)).*...
    sa)*2*pi/N;
end %getRS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fc,lambda] = getColForce(~,A,b,x0,jaco)

tol_rel = 0.0;
tol_abs = 0.0;
max_iter = 50;
[lambda, ~, ~, ~, ~, ~] = fischer_newton(A, b, x0, ...
                max_iter, tol_rel, tol_abs, 'perturbation', false);
            
fc = jaco'*lambda;
end % getColForce

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ns,totalnv,Xstart,Xend,totalPts] = preColCheck(o,X0,X1,walls,...
                                                upSampleFactor)
np = o.num_particles;                                            
nw = o.num_walls;
Np = o.points_per_particle;
Nw = o.points_per_wall;  

%Ns = ones(nv,1)*N*upSampleFactor;
Ns = ones(np,1)*Np;
totalnv = np;

% upsample positions
Xv0 = reshape(X0,Np,2*np);
% %Xv0 = interpft(Xv0,N,1);
% %Xv0 = interpft(Xv0,N*upSampleFactor,1);
Xv1 = reshape(X1,Np,2*np);
%Xv1 = interpft(Xv1,N*upSampleFactor,1);
%Xv1 = interpft(Xv1,N,1);

Xstart = Xv0(:);
Xend = Xv1(:);

if nw > 0
  Ns = [Ns;ones(nw,1)*Nw*upSampleFactor];
  totalnv = totalnv + nw;
  
  Xb = reshape(walls.X,Nw,2*nw);
  Xb = interpft(Xb,Nw*upSampleFactor,1);
  Xstart = [Xstart;Xb(:)];
  Xend = [Xend;Xb(:)];
end

totalPts = sum(Ns);

end % preColCheck

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vInf = bgFlow(~,X,options)
    
Np = size(X,1)/2;
np = size(X,2);
oc = curve;

[x,y] = oc.getXY(X);

if options.confined
    switch options.far_field
        case 'constant'
            vInf = [ones(Np,np);zeros(Np,np)];

        case 'circle'
            vInf = [-y(:,1); x(:,1)];
            
        case 'couette'

            vInf = zeros(2*Np,1);
            vInf = [vInf, [-options.couette_speed*y(:,2);...
                            options.couette_speed*x(:,2)]];
            
        case 'pipe'            
            vInf = [1 - y(:,1).^2; zeros(Np,1)];
            
        case 'bounded_shear'            
            if np == 1
                vInf = [-y(:,1); zeros(Np,1)];
            else if np == 2
                    vInf = [[-y(:,1);zeros(2*Np,1)],[-y(:,2);zeros(2*Np,1)]];
                end
            end
        otherwise
            vInf = [zeros(Np,np);zeros(Np,np)];
    end
else
    switch options.far_field
        case 'shear'
            vInf = [-3*X(end/2+1:end,:);zeros(Np,np)];

        case 'extensional'
            vInf = [X(1:end/2,:);-X(end/2+1:end,:)];

        case 'pipe'
            vInf = [1 - X(end/2+1:end,:).^2; zeros(Np,np)];
        otherwise
            vInf = [ones(Np,np);zeros(Np,np)];
    end

end
    
end % bgFlow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = build_matrix(o, geom, walls)
    
    np = o.num_particles;
    nw = o.num_walls;
    Np = o.points_per_particle;
    Nw = o.points_per_wall;
    
    Ntotal = 2*Np*np + 3*np;
    if nw > 0
        Ntotal = Ntotal + 2*Nw*nw +3*(nw-1);
    end
    
    I = eye(Ntotal);
    M = zeros(Ntotal);
    
    for k = 1:Ntotal
       M(:,k) = o.timeMatVec(I(:,k),geom,walls,true);
    end    
end % build_matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = build_preconditioner_matrix(o, geom, walls)

    [Np,np] = size(geom.X);
    [Nw,nw] = size(walls.X);
    
    Np = Np/2;
    Nw = Nw/2;
    
    Ntotal = 2*Np*np + 2*Nw*nw + 3*np + 3*(nw-1);
    
    M = zeros(Ntotal,Ntotal);
    
    % fibre-fibre
    for k = 1:np
        start_row = 1+(k-1)*2*Np;
        M(start_row:start_row+2*Np-1,start_row:start_row+2*Np-1) ...
                            =  -1/2*eye(2*Np)+o.Df(:,:,k);
    end
    
    % wall-wall
    for k = 1:nw
        start_row = 2*Np*np+1+(k-1)*2*Nw;
        M(start_row:start_row+2*Nw-1,start_row:start_row+2*Nw-1) ...
                            = -1/2*eye(2*Nw)+o.Dw(:,:,k)+o.N0w(:,:,k);
    end
    
    % u and omega
    for k = 1:np
        
        start_row = 2*Np*np+2*Nw*nw+2*(k-1)+1;
        start_col = 2*Np*(k-1)+1;
        
        M(start_row:start_row+1,start_col:start_col+2*Np-1) = ...
            [[-geom.sa(:,k)'*2*pi/Np, zeros(1,Np)];...
            [zeros(1,Np), -geom.sa(:,k)'*2*pi/Np]];
        
        start_row = 2*Np*np+2*Nw*nw+2*np+k;
        M(start_row,start_col:start_col+2*Np-1) = ...
            [(-geom.X(end/2+1:end,k)'+geom.center(2,k)).*geom.sa(:,k)'*2*pi/Np ...
            (geom.X(1:end/2,k)'-geom.center(1,k)).*geom.sa(:,k)'*2*pi/Np];
        
        start_row = 2*Np*(k-1)+1;
        start_col = 2*Np*np+2*Nw*nw+2*(k-1)+1;
        
        M(start_row:start_row+2*Np-1,start_col:start_col+1) = ...
            [[-ones(Np,1);zeros(Np,1)] ...
            [zeros(Np,1);-ones(Np,1)]];
        
        start_col = 2*Np*np + 2*Nw*nw+2*np+k;
        M(start_row:start_row+2*Np-1,start_col) = ...
            [geom.X(end/2+1:end,k)-geom.center(2,k); ...
            -geom.X(1:end/2,k)+geom.center(1,k)];
    end
    
    
    % stokeslets and rotlets
    oc = curve;
    sa = walls.sa;
    [x,y] = oc.getXY(walls.X);
    [cx,cy] = oc.getXY(walls.center);
    
    for k = 1:nw-1
        r = [x(:,k+1) - cx(k+1), y(:,k+1) - cy(k+1)];
        rho2 = (x(:,k+1) - cx(k+1)).^2 + (y(:,k+1) - cy(k+1)).^2;
        
        col_stokes1 = [-0.5*log(rho2) + r(:,1).*r(:,1)./rho2; ...
                        r(:,2).*r(:,1)./rho2]/(4*pi);
        col_stokes2 = [r(:,2).*r(:,1)./rho2; -0.5*log(rho2)...
                        + r(:,2).*r(:,2)./rho2]/(4*pi);
        col_rot = [r(:,2)./rho2; -r(:,1)./rho2];
        
        int_stokes = 2*pi*sa(:,k+1)'/Nw;
        int_rot = 2*pi*[(y(:,k+1).*sa(:,k+1))', -(x(:,k+1).*sa(:,k+1))']/Nw;
        
        start_row = 2*Np*np+2*Nw*nw+3*np+1;
        start_col = 2*Np*np+2*Nw*k+1;
        
        M(start_row:start_row+2,start_col:start_col+2*Nw-1) = ...
                    [int_stokes, zeros(1,Nw);...
                    zeros(1,Nw), int_stokes;...
                    int_rot];
        
        start_col = 2*Np*np+2*Nw*nw+3*np+(k-1)+1;
        M(start_row:start_row+2,start_col:start_col+2) = -2*pi*eye(3);
        
        start_row =  2*Np*np+2*Nw*k+1;
        start_col = 2*Np*np+2*Nw*nw+3*np+1 +3*(k-1);
        
        M(start_row:start_row+2*Nw-1,start_col:start_col+2) ...
                            = [col_stokes1, col_stokes2, col_rot];
    end
    
end % build_preconditioner_matrix

end % methods

end % classdef


