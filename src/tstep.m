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
Dp_old                  % Stokes double-layer potential for fiber-fiber interaction from previous time step
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
debug
explicit
matvecs

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
o.explicit = options.explicit;

o.resolve_collisions = options.resolve_collisions;
o.dt = prams.T/prams.number_steps;                  
o.gmres_tol = options.gmres_tol;             
o.gmres_max_it = 2*(prams.np*prams.Np + prams.nw*prams.Nw);

o.far_field = @(X) o.bgFlow(X,options); 
o.om = om;
o.tau0 = tau0;

o.num_particles = prams.np;
o.num_walls = prams.nw;
o.points_per_particle = prams.Np;
o.points_per_wall = prams.Nw;
o.debug = options.debug;
o.matvecs = 0;

np = o.num_particles;
nw = o.num_walls;
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
    o.Dp_old = zeros(size(o.Dp));
else
    o.Dp = [];
end
% CREATE UPSAMPLED MATRICES
% FIBRE-FIBRE
if ~isempty(geom.X)
    o.Dup = o.Dp;
else
    o.Dup = [];
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
            [-(geom.X(end/2+1:end,k)-geom.center(2,k)); ...
                (geom.X(1:end/2,k)-geom.center(1,k))];...
            [geom.sa(:,k)'*2*pi/Np/(2*pi) zeros(1,Np) 0 0 0];
            [zeros(1,Np) geom.sa(:,k)'*2*pi/np/(2*pi) 0 0 0];
            [(geom.X(end/2+1:end,k)'-geom.center(2,k)).*geom.sa(:,k)'*2*pi/Np ...
            -(geom.X(1:end/2,k)'-geom.center(1,k)).*geom.sa(:,k)'*2*pi/Np 0 0 0]/(2*pi)]);
    end
    
    if o.confined
        
        % WALL-WALL PRECONDITIONER
        o.precow.L = zeros(2*Nw + 3,2*Nw + 3,Nw);
        o.precow.U = zeros(2*Nw + 3,2*Nw + 3,Nw);
        
        oc = curve;
        sa = walls.sa;
        [x,y] = oc.getXY(walls.X);
        [cx,cy] = oc.getXY(walls.center);
        
        for k = 1:nw
            
            if k == 1 % first wall does not need Rotlets and Stokeslets
                [o.precow.L(1:2*Nw,1:2*Nw,k),o.precow.U(1:2*Nw,1:2*Nw,k)] =...
                    lu(-1/2*eye(2*Nw)+o.Dw(:,:,k)+o.N0w(:,:,k));
                
            else
                r = [x(:,k) - cx(k), y(:,k) - cy(k)];
                rho2 = (x(:,k) - cx(k)).^2 + (y(:,k) - cy(k)).^2;
                
                col_stokes1 = [-0.5*log(rho2) + r(:,1).*r(:,1)./rho2; ...
                                r(:,2).*r(:,1)./rho2]/(4*pi);
                col_stokes2 = [r(:,2).*r(:,1)./rho2; ...
                                -0.5*log(rho2) + r(:,2).*r(:,2)./rho2]/(4*pi);
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
                uP_m1,omega_m1,first_step,forceP_old, torqueP_old, etaP_old, ...
                etaW_old)
% timeStep(geom, walls, xc, tau) takes a single time step and returns the
% angles and velocities of the particles at the next time step as well as
% the density functions on the particles and walls, the translational and
% rotational velocities of the particles, and the net force and torque on
% all particles and walls

np = o.num_particles;
Np = o.points_per_particle;
nw = o.num_walls;

% ASSEMBLE RHS WITH NO FORCE AND TORQUE ON PARTICLES
forceP = zeros(2,np);
torqueP = zeros(1,np);
forceW = zeros(2,nw-1);
torqueW = zeros(1,nw-1);

explicit_time_step = o.explicit; % && ~first_step;

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
dtau = tau - o.tau0;
o.tau0 = tau;

o.Dp_old = o.Dp;    

% ROTATE EXISTING PRECONDITIONER
for i = 1:np

    R = spdiags([sin(dtau(i))*ones(2*Np,1), cos(dtau(i))*ones(2*Np,1)...
                -sin(dtau(i))*ones(2*Np,1)], [-Np, 0, Np], zeros(2*Np, 2*Np));
          
    o.Dp(:,:,i) = R*o.Dp_old(:,:,i)*R';
    o.Dup(:,:,i) = o.Dp(:,:,i);
    
    if o.use_precond
        
        [o.precop.L(:,:,i),o.precop.U(:,:,i)] =...
            lu([-1/2*eye(2*Np)+o.Dp(:,:,i) ...
            [-ones(Np,1);zeros(Np,1)] ...
            [zeros(Np,1);-ones(Np,1)] ...
            [-(geom.X(end/2+1:end,i)-geom.center(2,i));...
                    geom.X(1:end/2,i)-geom.center(1,i)];...
            [geom.sa(:,i)'*2*pi/Np/(2*pi) zeros(1,Np) 0 0 0];
            [zeros(1,Np) geom.sa(:,i)'*2*pi/Np/(2*pi) 0 0 0];
            [(geom.X(end/2+1:end,i)'-geom.center(2,i)).*geom.sa(:,i)'*2*pi/Np ...
            -(geom.X(1:end/2,i)'-geom.center(1,i)).*geom.sa(:,i)'*2*pi/Np 0 0 0]/(2*pi)]);
    end    
end

if explicit_time_step
    rhs = o.assembleRHS_explicit(geom, walls, forceP, torqueP, ...
                    etaP_old, 1:np, false);
else
    rhs = o.assembleRHS(geom, walls, forceP, torqueP, forceW, torqueW, 1:np, false);
end

% SOLVE SYSTEM USING GMRES TO GET CANDIDATE TIME STEP
if o.use_precond
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom,walls,1,false),rhs,[],...
      o.gmres_tol,o.gmres_max_it,@o.preconditionerBD,[]);
else
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom,walls,1,false),rhs,[],...
      o.gmres_tol, o.gmres_max_it, [], []);
end

iter = I(2);

% UNPACK SOLUTION
[etaP, etaW, uP, omega, forceW, torqueW] = o.unpackSolution(Xn,geom,walls,true);
omega = -omega;

% CREATE CANDIDATE CONFIGURATION 
if first_step || o.tstep_order == 1 % FORWARD EULER
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
    geomOld = geom;
    
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
    solution.forceP = zeros(2,np);
    solution.torqueP = zeros(1,np);
    solution.forceP_old = forceP_old;
    solution.torqueP_old = torqueP_old;
    solution.etaP_old = etaP_old;
    solution.etaW_old = etaW_old;
    
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

end % timeStep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tx = timeMatVec(o, Xn, geom, walls, bounding_wall, preprocess)
% Tx = timeMatVec(Xn,geom) does a matvec for GMRES 
if o.profile
    tMatvec = tic;
end

o.matvecs = o.matvecs + 1;

np = geom.n;
Np = geom.N;

if ~isempty(walls)
    nw = walls.n;
    Nw = walls.N;
else
    nw = 0;
    Nw = 0;
end

pot_walls = poten(Nw, o.om);
pot_particles = poten(Np, o.om);

if np ~= o.num_particles
    Dp = pot_particles.stokesDLmatrix(geom);
else
    Dp = o.Dp;
end

if nw ~= o.num_walls && nw > 0
    Dw = pot_particles.stokesDLmatrix(walls);
else
    Dw = o.Dw;
end

% Xsou = geom.X; 
% Nup = Np*ceil(sqrt(Np));
% 
% Xup = [interpft(Xsou(1:Np,:),Nup);...
%    interpft(Xsou(Np+1:2*Np,:),Nup)];
% 
% geomUp = capsules([],Xup);
% Dup =  pot_particles.stokesDLmatrix(geomUp);
    
% PREALLOCATE OUTPUT VECTORS
% velParticles : velocity of particles
% velWalls     : velocity of particles
% forceP       : net force on particles, strength of Stokeslets
% torqueP      : net torque on particles, strengh of rotlets
velParticles = zeros(2*Np,np);
if nw > 0
    velWalls = zeros(2*Nw,nw);
else
    velWalls = [];
end

% UNPACK Xn
% etaP    : denisty function on particles
% etaW    : density function on walls
% uT      : translational velocity of particles
% omega   : angular velocity of particles
% forceW  : net force on walls, strength of Stokeslets
% torqueW : net torque on walls, strength of rotlets
[etaP, etaW, uP, omega, forceW, torqueW] = unpackSolution(o, Xn, geom, walls, bounding_wall ~= 0);

% CALCULATE VELOCITIES ON PARTICLES AND WALLS, NET FORCE AND TORQUE ON
% WALLS

% pp_dlp : velocity on particles due to density funtion on other particles
% wp_dlp : velocity on walls due to density function on all particles
% pw_dlp : velocity on particles due to density function on all walls
% ww_dlp : velocity on walls due to denisity function on other walls
% p_sr   : velocity on particles due to rotlets and Stokeslets in walls
% w_sr   : velocity on walls due to rotlets and Stokeslets in walls

% ADD JUMP IN DLP
velParticles = velParticles - 1/2*etaP;

if o.confined && nw > 0
   velWalls = velWalls - 1/2*etaW;
end 

% ADD SELF CONTRIBUTION
velParticles = velParticles + ...
    pot_particles.exactStokesDLdiag(geom, Dp, etaP);

if o.confined && nw > 0
    velWalls = velWalls + pot_walls.exactStokesDLdiag(walls, o.Dw, etaW);
end

% START OF SOURCE == PARTICLES
% START OF TARGET == PARTICLES

if ~o.explicit
    if o.fmm
        kernel = @pot_particles.exactStokesDLfmm;
    else
        kernel = @pot_particles.exactStokesDL;
    end
    
    kernelDirect = @pot_particles.exactStokesDL;
    
    if o.near_singular
        near_structff = geom.getZone([],1);
        
        DLP = @(X) pot_particles.exactStokesDLdiag(geom,Dp,X) - 1/2*X;
        pp_dlp = pot_particles.nearSingInt(geom, etaP, DLP, Dp, ...
            near_structff ,kernel, kernelDirect, geom, true, false);
    else
        pp_dlp = kernel(geom, etaP, Dp);
    end
else
    pp_dlp = zeros(2*Np,np);
end
% END OF TARGET == PARTICLES

% START OF TARGET == WALLS 
if o.confined && np > 0  && nw > 0
   
   if o.fmm
       kernel = @pot_walls.exactStokesDLfmm;       
   else
       kernel = @pot_walls.exactStokesDL;
   end
   
   kernelDirect = @pot_walls.exactStokesDL;
   
   if o.near_singular
       DLP = @(X) pot_walls.exactStokesDLdiag(geom, Dp, X) - 1/2*X;
       wp_dlp = pot_walls.nearSingInt(geom, etaP, DLP, Dp, ...
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
if o.confined && np > 0  && nw > 0
   
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
if o.confined && nw > 0
    
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

if o.confined && nw > 0
    % START SOURCE == ROTLETS AND STOKESLETS
    % START TARGETS == PARTICLES
    p_sr = 0;
    
    if bounding_wall ~= 0
        start_index = 2;
    else
        start_index = 1;
    end
    
    for k = start_index:nw % loop over all walls, except outer wall
            p_sr = p_sr + ...
                o.computeRotletStokesletVelocity(geom.X, walls.center(:,k), ...
                forceW(:,k), torqueW(k));
    end
    
    % END TARGETS == PARTICLES
    
    % START TARGETS == WALLS
    w_sr = zeros(2*Nw,nw);
    for k = start_index:nw % loop over all walls, except outer wall
        if k ~= bounding_wall
            w_sr = w_sr +...
                o.computeRotletStokesletVelocity(walls.X, walls.center(:,k),...
                                    forceW(:,k), torqueW(k));
        end
    end   
    
    % END TARGETS == WALLS
    % END SOURCE == ROTLETS
    
    % COMPUTE NET FORCE AND TORQUE ON WALLS
    [f_w, t_w] = o.computeNetForceTorque(etaW, walls);
    
    if preprocess
        forceW = f_w; 
        torqueW = t_w;
    else
        forceW = f_w - forceW; 
        torqueW = t_w - torqueW;
    end
else
    p_sr = 0;
    w_sr = 0;
    forceW = [];
    torqueW = [];
end

% EVALUATE TOTAL VELOCITY ON PARTICLES

% ADD CONTRIBUTIONS FROM OTHER BODIES
velParticles = velParticles + p_sr + pp_dlp + pw_dlp;

% SUBTRACT VELOCITY ON SURFACE
for k = 1:np
  velParticles(1:Np,k) = velParticles(1:Np,k) - uP(1,k);
  velParticles(Np+1:end,k) = velParticles(Np+1:end,k) - uP(2,k);
end

for k = 1:np 
  velParticles(1:Np,k) = velParticles(1:Np,k) ...
                - (geom.X(Np+1:end,k) - geom.center(2,k))*omega(k);
  velParticles(Np+1:end,k) = velParticles(Np+1:end,k)...
                + (geom.X(1:Np,k) - geom.center(1,k))*omega(k); 
end

% EVALUATE VELOCITY ON WALLS
if o.confined && nw > 0 && bounding_wall ~= 0
    velWalls = velWalls + ww_dlp + wp_dlp + w_sr;
    velWalls(:,1) = velWalls(:,1)...
                    + pot_walls.exactStokesN0diag(walls, o.N0w, etaW(:,1));
end

% EVALUTATE FORCES ON PARTICLES
[forceP, torqueP] = o.computeNetForceTorque(etaP, geom);

% CONSTRUCT OUTPUT VECTOR
if nw > 0
    
    if bounding_wall == 0
        Tx = [velParticles(:); velWalls(:); forceP(:); torqueP(:); ...
                        forceW(:); torqueW(:)];
    else
        Tx = [velParticles(:); velWalls(:); forceP(:); torqueP(:); ...
                        forceW(:,2:end); torqueW(2:end)];
    end
else
    Tx = [velParticles(:); forceP(:); torqueP(:)];
end

if o.profile
    o.om.writeMessage(['Matvec assembly completed in ', ....
                            num2str(toc(tMatvec)), ' seconds']);
end

end % timeMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhs = assembleRHS(o, geom, walls, forceP, torqueP, forceW, torqueW, ...
                            bodies, preprocess)
% Assemble right hand side
np = geom.n;
Np = geom.N;

if ~isempty(walls)
    nw = walls.n;
    Nw = walls.N;
else
    nw = 0;
    Nw = 0;
end

if ~preprocess
    if o.confined
        ff = o.far_field(walls.X);
    else
        ff = -o.far_field(geom.X); 
    end
end

if o.confined
    if ~preprocess
        rhs = [zeros(2*Np*np,1); ff(:); forceP(:); ...
                    torqueP'; forceW(:); torqueW'];
    else
        rhs = [zeros(2*Np*np,1); zeros(2*Nw*nw,1); forceP(:); torqueP'; ...
                            forceW(:); torqueW'];
    end
else
    if ~preprocess
        rhs = [ff(:); forceP(:); torqueP(:)];
    else
        rhs = [zeros(2*Np*np,1); forceP(:); torqueP(:)];
    end
end

if o.resolve_collisions
 
    % ADD IN CONTRIBUTIONS TO VELOCITIES FROM PARTICLE ROTLETS AND STOKESLETS
    for k = bodies
	
        % CONTRIBUTION TO PARTICLE VELOCITY
        v = o.computeRotletStokesletVelocity(geom.X, geom.center(:,k),...
            	forceP(:,k),torqueP(k));

	    rhs(1:2*Np*np) = rhs(1:2*Np*np) - v(:);

        if o.confined && ~isempty(walls)
            % CONTRIBUTION TO WALL VELOCITY FROM PARTICLE SINGULARITIES
            v = o.computeRotletStokesletVelocity(walls.X, geom.center(:,k),...
                                forceP(:,k),torqueP(k));

            rhs(2*Np*np+1:2*Np*np+2*Nw*nw) = ...
                                rhs(2*Np*np+1:2*Np*np+2*Nw*nw) - v(:);
            
        end
    end
    
    if o.confined && ~isempty(walls)
        
        % CONTRIBUTION TO WALL AND PARTICLE VELOCITY FROM WALL SINGULARITIES        
        for k = 1:length(torqueW)
            
            v = o.computeRotletStokesletVelocity(geom.X, walls.center(:,k),...
                forceW(:,k),torqueW(k));
            
            rhs(1:2*Np*np) = rhs(1:2*Np*np) - v(:);
            
            v = o.computeRotletStokesletVelocity(walls.X, walls.center(:,k),...
                forceW(:,k),torqueW(k));
            
            rhs(2*Np*np+1:2*Np*np+2*Nw*nw) = ...
                rhs(2*Np*np+1:2*Np*np+2*Nw*nw) - v(:);
        end
    end
    
end

end % assembleRHS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhs = assembleRHS_explicit(o, geom, walls, forceP, torqueP, ...
                    etaP_old, bodies, preprocess)
% Assemble right hand side

np = geom.n;
nw = o.num_walls;
Np = geom.N;
Nw = o.points_per_wall;

if ~preprocess
    if o.confined
        ff = o.far_field(walls.X);
    else
        ff = -o.far_field(geom.X); 
    end
end

if o.confined
    if ~preprocess
        rhs = [zeros(2*Np*np,1); ff(:); forceP(:); ...
                    torqueP'; zeros(3*(nw-1),1)];
    else
        rhs = [zeros(2*Np*np,1); forceP(:); torqueP'];
    end
else
    if ~preprocess
        rhs = [ff(:); forceP(:); torqueP(:)];
    else
        rhs = [zeros(2*Np*np,1); forceP(:); torqueP(:)];
    end
end

% APPLY VELCOITY FROM OLD DENSITY FUNCTION TO ALL OTHER BODIES
if ~preprocess
    pot_particles = o.potp;
    if o.fmm
        kernel = @pot_particles.exactStokesDLfmm;
    else
        kernel = @pot_particles.exactStokesDL;
    end
    
    kernelDirect = @pot_particles.exactStokesDL;
    
    if o.near_singular
        DLP = @(X) pot_particles.exactStokesDLdiag(geom,o.Dp_old,X) - 1/2*X;
        pp_dlp = pot_particles.nearSingInt(geom, etaP_old, DLP, o.Dp_old, ...
            o.near_structff, kernel, kernelDirect, geom, true, false);
    else
        pp_dlp = kernel(geom, etaP_old, o.Dp_old);
    end
    
    rhs(1:2*np*Np) = rhs(1:2*np*Np) - pp_dlp(:);
end

if o.resolve_collisions
    
    for k = bodies
        
        % CONTRIBUTION TO PARTICLE VELOCITY
        v = o.computeRotletStokesletVelocity(geom.X, geom.center(:,k),...
            forceP(:,k),torqueP(k));
        
        % SUBTRACT VELOCITY DUE TO CURRENT FORCE ON SELF
        rhs(2*Np*(k-1)+1:2*Np*k) = rhs(2*Np*(k-1)+1:2*Np*k) - v(:,k);
        
        if ~preprocess
            
            % SUBTRACT OFF VELOCITY FROM ROTLETS INDUCED BY etaP_old
            [force_tmp, torque_tmp] = o.computeNetForceTorque(etaP_old, geom);

            v = o.computeRotletStokesletVelocity(geom.X, geom.center(:,k),...
                force_tmp(:,k),torque_tmp(k));
            
            % intra-particle interactions have already been considered
            v(:,k) = 0;
            
            rhs(1:2*Np*np) = rhs(1:2*Np*np) - v(:);
            
            if o.confined
                % CONTRIBUTION TO WALL VELOCITY
                v = o.computeRotletStokesletVelocity(walls.X, geom.center(:,k),...
                    forceP(:,k),torqueP(k));
                rhs(2*Np*np+1:2*Np*np+2*Nw*nw) = ...
                    rhs(2*Np*np+1:2*Np*np+2*Nw*nw) - v(:);
            end
        end
    end
end

end % assembleRHS_explicit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [etaP, etaW, u, omega, forceW, torqueW] = unpackSolution(...
                            o, Xn, geom, walls, bounding_wall)
% Unpack solution vector that comes from GMRES call

Np = o.points_per_particle;
Nw = o.points_per_wall;

np = geom.n;

if ~isempty(walls)
    nw = walls.n;
else
    nw = 0;
end

% REORGANIZE COLUMN VECTOR INTO MATRIX
% EXTRACT DENSITY FUNCITONS ON PARTICLES AND WALLS
% each column of etaP corresponds to the density function of a rigid body
etaP = zeros(2*Np,np);
for k = 1:np
  etaP(:,k) = Xn((k-1)*2*Np+1:k*2*Np);
end

% each column of etaW corresponds to the density function of a solid wall
if o.confined && nw > 0
    etaW = zeros(2*Nw, nw);
    for k = 1:nw
       etaW(:,k) = Xn(2*Np*np+1+(k-1)*2*Nw:2*Np*np+k*2*Nw);
    end
else
    etaW = 0;
end

% EXTRACT TRANSLATIONAL AND ROTATIONAL VELOCITIES
u = zeros(2,np);
omega = zeros(1,np);
if nw > 0
    for k = 1:np
        u(:,k) = Xn(2*Np*np+2*Nw*nw+(k-1)*2+1:2*Np*np+2*Nw*nw+k*2);
        
        omega(k) = Xn(2*Np*np+2*Nw*nw+2*np+k);
    end
else
    for k = 1:np
        u(:,k) = Xn(2*Np*np+(k-1)*2+1:2*Np*np+k*2);
        
        omega(k) = Xn(2*Np*np+2*np+k);
    end
end

if o.confined && nw > 0
    
    forceW = zeros(2,nw);
    if bounding_wall
        end_index = nw - 1;
    else
        end_index = nw;
    end

    for k = 1:end_index
       forceW(:,k) = ...
           Xn(2*Np*np+2*Nw*nw+3*np+1+2*(k-1): 2*Np*np+2*Nw*nw+3*np+2*k);
    end

    torqueW = zeros(1,nw);
    for k = 1:end_index
       torqueW(k) = Xn(2*Np*np+2*Nw*nw+3*np+2*end_index+k);
    end
    
    % reorganize to put bounding wall first if needed
    if bounding_wall
       forceW = [forceW(:,end), forceW(:,1:end-1)];
       torqueW = [torqueW(end), torqueW(1:end-1)];
    end
else
    forceW = 0;
    torqueW = 0;
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
function vel = computeRotletStokesletVelocity(~,X,center,stokeslet,rotlet)
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

RotTerm = ((y-cy)./rho2)*rotlet;
velx = (LogTerm + rorTerm)/(4*pi) + RotTerm/(4*pi);
%velx = 1/4/pi*(LogTerm + rorTerm + RotTerm);
% x component of velocity due to the stokeslet and rotlet

LogTerm = -0.5*log(rho2)*stokeslet(2);
rorTerm = 1./rho2.*((y-cy).*(x-cx)*stokeslet(1) + ...
                 (y-cy).*(y-cy)*stokeslet(2));

RotTerm = -((x-cx)./rho2)*rotlet;
vely = (LogTerm + rorTerm)/(4*pi) + RotTerm/(4*pi);
%vely = 1/4/pi*(LogTerm + rorTerm + RotTerm);
% y component of velocity due to the stokeslet and rotlet

vel = [velx;vely];

end % computeRotletStokesletVelocity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [force,torque] = computeNetForceTorque(~, eta, geom)
% z = letsIntegrals(stokeslet,rotlet,etaM,walls) integrates the density
% function to enforce constraints on stokeslets and rotlets

N = geom.N;
n = geom.n;

force = zeros(2,n);
torque = zeros(1,n);

for k = 1:n

  force(1,k) =  sum(eta(1:N,k).*geom.sa(:,k))*2*pi/N;
  force(2,k) =  sum(eta(N+1:end,k).*geom.sa(:,k))*2*pi/N;

  torque(k) = sum(((geom.X(N+1:end,k) - geom.center(2,k)).*eta(1:N,k) - ...
    (geom.X(1:N,k) - geom.center(1,k)).*eta(N+1:end,k)).*geom.sa(:,k))*2*pi/N;

  force(:,k) = force(:,k)/(2*pi);
  torque(k) = torque(k)/(2*pi);  
end

end % computeNetForceTorque

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solution = resolveCollisions(o, geomOld, geomProv, walls, ...
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
%   [~,length] = oc.geomProp(walls.X);
%   walllength = max(length/Nw);
  %wallupSamp = ceil(walllength/cellSize);
  wallupSamp = 1;
else
  wallupSamp = 1;
end

upSampleFactor = 1;
nexten = 0;
c_tol = 1e-12;
minSep = o.minimum_separation;
maxIter = 1000;

% load explicit_dv1
% Xs = Xstart;
% Xf = Xend;
%     
% CHECK FOR COLLISION
[Ns,totalnp,Xstart,Xend,totalPts] = o.preColCheck(geomOld.X,geomProv.X,...
                                        walls,wallupSamp);
                                     
[vgrad, iv, ids, vols] = getCollision(Ns, totalnp, Xstart, Xend, minSep, ...
        maxIter, totalPts, c_tol, np, nw, Np*upSampleFactor, ...
        Nw*wallupSamp, nexten, max(cellSize,minSep));
    
o.om.writeMessage(['ivolume: ' num2str(iv)],'%s\n');

colCount = 0;

forceP = zeros(2, np);
torqueP = zeros(1,np);

alpha0 = 0.75;
alpha = alpha0;
alpha_max = 10;
alpha_min = 0.25;

% RESOLVE COLLISION
while(iv < 0)

    colCount = colCount + 1;  
    o.om.writeMessage(['Collision resolving iteration: ', num2str(colCount)]);
    
    %vgrad = vgrad(1:2*Np*np);
% 
%     ids = ids(1:2*Np*np);
%     vols = vols(1:2*Np*np);

                            
    [A,~,ivs,~,jacoSmooth, vtoiv] = o.preprocessRigid(vgrad,ids,vols,geomOld,...
                    walls);
    
    % CALCULATE FORCE AND TORQUES ON PARTICLES
    lambda = -A\(ivs/o.dt); 
    
    % if all lambda values are negative, adjust vgrad and try again
    if max(lambda) < 0
        vgrad(1:2*Np*np) = o.adjustNormal(vgrad(1:2*Np*np),Np,np,geomOld,...
            edgelength);
        [A,~,ivs,~,jacoSmooth, vtoiv] = o.preprocessRigid(vgrad,ids,vols,geomOld,...
            walls);
        
        lambda = -A\(ivs/o.dt);
    end
    
    fc = jacoSmooth'*lambda;
    
    % if at least one of the lambda values is still negative use Fisher-Newton to
    % solve for best positive lambdas
    if min(lambda)< 0
        [fc, ~] = o.getColForce(A,ivs/o.dt,-lambda,jacoSmooth);
    end    
    
    if colCount > 2
        alpha = max(min(o.compute_alpha(alpha, alpha0, iv), alpha_max), alpha_min); 
        disp(['alpha = ', num2str(alpha)]);              
    end
    
    fc = alpha*fc;
    
    %fc_tmp = o.f_smooth(full(fc),Np,np);
    %fc_tot = fc_tot + fc;
    
    %fc_tot_tmp = reshape(fc_tot, 2*Np, np);
    fc_tmp_particles = reshape(fc(1:2*Np*np), 2*Np, np);
    [forceP_tmp, torqueP_tmp] = o.computeNetForceTorque(fc_tmp_particles,geomOld);
    
    if o.confined
        fc_tmp_walls = reshape(fc(2*Np*np + 1:end), 2*Nw, nw);
        [forceW_tmp, torqueW_tmp] = o.computeNetForceTorque(fc_tmp_walls,walls);
        
        forceW = forceW_tmp(:,2:end);
        torqueW = torqueW_tmp(2:end);
    end
    
    forceP = forceP + o.balance_force(geomOld,forceP_tmp, vtoiv);    
    torqueP = torqueP  + torqueP_tmp;
    
    % ASSEMBLE NEW RHS WITH CONTACT FORCES
    if o.explicit
        rhs = o.assembleRHS_explicit(geomOld, walls, forceP, torqueP, ...
            solution.etaP_old, 1:np, false);
    else
        rhs = o.assembleRHS(geomOld, walls, forceP, torqueP, forceW, torqueW, 1:np, false);
    end
    
    % SOLVE SYSTEM
    if o.use_precond
        Xn = gmres(@(X) o.timeMatVec(X,geomOld,walls,1,false),rhs,[],o.gmres_tol,...
            o.gmres_max_it,@o.preconditionerBD,[]);
    else
        Xn = gmres(@(X) o.timeMatVec(X,geomOld,walls,1,false),rhs,[],o.gmres_tol,...
            o.gmres_max_it);
    end 

    % UNPACK SOLUTION
    [etaP, etaW, uP, omega, forceW, torqueW] = o.unpackSolution(Xn, geomOld, walls, true);
 
    if o.debug
        close all
        hold on
        v = zeros(2*Np,np);

        u_max_old = norm(max(abs(solution.uP), [], 2));
        f_max = norm(max(abs(forceP), [], 2));
        
        for k = 1:np
            plot(geomOld.X(1:end/2,k),geomOld.X(end/2+1:end,k),'-b');
            quiver(geomOld.X(1:end/2,k),geomOld.X(end/2+1:end,k), ...
                vgrad((k-1)*2*Np + 1:(k-1)*2*Np+Np),...
                vgrad((k-1)*2*Np+Np+1:(k-1)*2*Np+2*Np),'r');
            quiver(geomOld.center(1,k),geomOld.center(2,k),uP(1,k)/u_max_old,uP(2,k)/u_max_old,1,'g');
            quiver(geomOld.center(1,k),geomOld.center(2,k),solution.uP(1,k)/u_max_old,solution.uP(2,k)/u_max_old,1,'c');
            quiver(geomOld.center(1,k),geomOld.center(2,k),forceP(1,k)/f_max,forceP(2,k)/f_max,1,'m');        

            v = v + o.computeRotletStokesletVelocity(geomOld.X, geomOld.center(:,k),...
                forceP(:,k),torqueP(k));
        end
        
        quiver(geomOld.X(1:end/2,:), geomOld.X(end/2+1:end,:), v(1:end/2,:), v(end/2+1:end,:));
        
        axis equal
        pause
     end
    
    omega = -omega;
    % UPDATE PARTICLE POSTIONS AND ANGLES 
    if first_step || o.tstep_order == 1% FORWARD EULER
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

    [vgrad, iv(end+1), ids, vols] = getCollision(Ns, totalnp, Xstart, Xend, minSep,... 
        maxIter, totalPts, c_tol, np, nw, Np*upSampleFactor, Nw*wallupSamp,...
        nexten, max(cellSize,minSep));

    solution.xc_new = xc_new;
    solution.tau_new = tau_new;
    solution.etaP = etaP;
    solution.etaW = etaW;
    solution.uP = uP;
    solution.omega = omega;
    solution.forceW = forceW;
    solution.torqueW = torqueW;
    [solution.forceP, solution.torqueP] = o.computeNetForceTorque(solution.etaP,geomOld);
    o.om.writeMessage(['ivolume: ' num2str(iv(end))],'%s\n');
end


end % resolveCollision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,jaco,ivs,listnv,jacoSmooth, vtoiv] = preprocessRigid(o,vgrad,ids,...
                    vols,geom,walls)
                
Np = o.points_per_particle;
np = o.num_particles;
nw = o.num_walls;
Nw = o.points_per_wall;

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

if o.confined
    for i = 2*Np*np + 1:length(ids)
        if (ids(i) ~= 0)
            k = ceil((i - 2*Np*np)/(2*Nw)) + np; % number walls after particles
            listnv = [listnv;k]; % keep running track of all walls
            jI = [ids(i);jI];    % keep track of all interference volumes
            jJ = [i;jJ];         % keep track of all nodes
            jV = [vgrad(i);jV];  % keep track of all volume gradients
            vtoiv(k,ids(i)) = 1;
            ivs(ids(i)) = vols(i);
        end
    end
end

listnv = unique(listnv);
ivs(ivs>-1e-12) = -1e-12;
% sparse matrix of size nivs by total number of points on particles
jaco = sparse(jI,jJ,jV,nivs,2*Np*np+2*Nw*nw); 
jacoSmooth = jaco*1;

% LOOP OVER INTERSECTION VOLUMES
for i = 1:nivs
    
    % FIND BODIES IN THIS VOLUME
    S = find(vtoiv(:,i)~=0); 
    
    fp = jaco(i,1:2*Np*np)';    
    fp = reshape(fp, 2*Np, np);

    [forceP,torqueP] = o.computeNetForceTorque(fp,geom);
    forceP = forceP(:,S(S<np)');
    torqueP = torqueP(S(S<np)');
    
    if o.confined
        
        fw = jaco(i,2*Np*np+1:end);
        fw = reshape(fw, 2*Nw, nw);
        [forceW,torqueW] = o.computeNetForceTorque(fw,walls);
        forceW = forceW(:,S(S>np) - np);
        torqueW = torqueW(S(S>np) - np);
    else
        forceW = [];
        torqueW = [];
    end
    
    geomTmp = capsules(o.prams, geom.X(:,S(S<np)'));

    % one of the "particles" in volume is a solid wall
    if max(S) > np
        wallsTmp = capsules(o.prams, walls.X(:, S(S > np) - np));
    else
        wallsTmp = [];
    end
 
    if max(S == np+1) % collision with outer wall
        bounding_wall = 1;
    else
        bounding_wall = 0;
    end
    
    % ASSEMBLE NEW RHS WITH CONTACT FORCES
    if o.explicit 
        rhs = o.assembleRHS_explicit(geomTmp, wallsTmp, forceP, torqueP, ...
                [], 1:length(S(S<np)), true);
    else
        rhs = o.assembleRHS(geomTmp, wallsTmp, forceP, torqueP, forceW, torqueW,...
                        1:length(S(S<np)), true);  
    end
    
    % SOLVE SYSTEM WITH FAR FIELD NEGLECTED
    
    Xn = gmres(@(X) o.timeMatVec(X,geomTmp,wallsTmp,bounding_wall,true),rhs,[],o.gmres_tol,4*Np);
    
    % COMPUTE VELOCITY OF EACH POINT ON ALL RIGID PARTICLES
    [~, ~, uP, omega, ~, ~] = o.unpackSolution(Xn, geomTmp, wallsTmp, bounding_wall == 1);

    for k = 1:geomTmp.n
        b = [uP(1,k)*ones(Np,1); uP(2,k)*ones(Np,1)] + ... %translational
            omega(k)*[geomTmp.X(Np+1:end,k) - geomTmp.center(2,k); ...% + rotational
            -(geomTmp.X(1:Np,k) - geomTmp.center(1,k))];
        
        j = S(k);
        SS = find(vtoiv(j,:)~=0);
        for l = 1:numel(SS)
            A(SS(l),i) = A(SS(l),i) + dot(jaco(SS(l),1+(j-1)*2*Np:2*Np+(j-1)*2*Np),b);
        end
    end

    if o.debug
        close all
        hold on
        for k = 1:geomTmp.n
            plot(geom.X(1:end/2,S(k)),geom.X(end/2+1:end,S(k)),'-bo');
            quiver(geom.X(1:end/2,S(k)),geom.X(end/2+1:end,S(k)), ...
                vgrad((S(k)-1)*2*Np + 1:(S(k)-1)*2*Np+Np),...
                vgrad((S(k)-1)*2*Np+Np+1:(S(k)-1)*2*Np+2*Np),'r');
            quiver(geomTmp.center(1,k),geomTmp.center(2,k),uP(1,k),uP(2,k),1000,'g');
            quiver(geomTmp.center(1,k),geomTmp.center(2,k),forceP(1,k),forceP(2,k),1000,'m');
            
            b = [uP(1,k)*ones(Np,1); uP(2,k)*ones(Np,1)] + ... %translational
                omega(k)*[geomTmp.X(Np+1:end,k) - geomTmp.center(2,k); ...% + rotational
                 -(geomTmp.X(1:Np,k) - geomTmp.center(1,k))];
            
           quiver(geomTmp.X(1:end/2,k),geomTmp.X(end/2+1:end,k), b(1:end/2),b(end/2+1:end),'k');
        end
        
        axis equal
        pause
    end
end

end % preprocessRigid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha = compute_alpha(~, alpha, alpha0, vols)

r = log(abs(vols(end)/vols(end-1)))/log(abs(vols(end-1)/vols(end-2)));
    
disp(['Convergence rate is ', num2str(r)]);
    
if r > 0
     alpha = alpha*4/(0.5*r + 1)^2;
else
    alpha = alpha0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = balance_force(o, geom, F, vtoiv)
    
%n_vols = size(vtoiv,2);

if o.debug
   
   disp(F);
   disp(vtoiv);
   close all
   hold on
   
   plot(geom.X(1:end/2,:),geom.X(end/2+1:end,:), 'b');
   quiver(geom.center(1,:),geom.center(2,:),F(1,:),F(2,:), 'r');
   axis equal
   
   pause
end

C = o.determine_clusters(vtoiv);

for ic = 1:length(C)
   
    % identify which particle is in the most volumes (i.e. has the most
    % neighbors)
    
    [~, max_ind] = max(sum(vtoiv(C{ic},:)'));
    target = C{ic}(max_ind);
    % total force on this particle will be negative sum of forces on all
    % other particles in cluser
    
    if length(C{ic}) > 1
        total_force = [0;0];
        for k = C{ic}'
            if k ~= target
                total_force = total_force + F(:,k);
            end
        end
    
        F(:,target) = -total_force;
    end
end

end % balance_force

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C, volumes] = add_particle_to_cluster(o, C, S, k, volumes)
    
n_vols = length(S);

for i = 1:n_vols
    if ~max(i == volumes)
        if max(S{i} == k)
            volumes(end+1) = i;
            C = [C; S{i}];

            for k1 = S{i}'
                [C, volumes] = o.add_particle_to_cluster(C, S, k1, volumes);
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = determine_clusters(o, vtoiv)
    
n_vols = size(vtoiv,2);
S = cell(n_vols,1);
C = cell(0);

for v = 1:n_vols
 
    % identify particles in this volume
    S{v} = find(vtoiv(:,v)~=0);    
    
end

% concatenate volumes
already_considered = 0;
% loop over volumes
for i = 1:n_vols
    
    % if this volume isn't already part of a cluster continue
    if ~max(i == already_considered)
        
        for k = S{i}'
            [C{end+1}, already_considered] = o.add_particle_to_cluster([], S, k, already_considered);
        end
    end
end

C = C(~cellfun('isempty', C));

for i = 1:length(C)
   C{i} = unique(C{i}); 
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vgrad] = adjustNormal(~,vgrad,N,n,geom,edgelength)
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
        
        %disp('Using surface normal for volume gradient');
        bnnorm = edgelength(k);
        vgrad_mag = sqrt(vgrad(i,k)^2 + vgrad(N+i,k)^2);
        
        tangent = geom.xt;
        oc = curve;
        [tanx,tany] = oc.getXY(tangent);
        nx = tany;
        ny = -tanx;
        normal = [nx;ny];
        
        normal_mag = sqrt(nx(i,k)^2 + ny(i,k)^2);
        
        vgrad(i,k) = vgrad_mag*normal(i,k)/normal_mag;
        vgrad(N+i,k) = vgrad_mag*normal(i+N,k)/normal_mag;
        
%         vgrad(i,k) = normal(i,k)*bnnorm;
%         vgrad(N+i,k) = normal(i+N,k)*bnnorm;
        
    end
end

vgrad = vgrad(:);
end % adjustNormal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fc,lambda] = getColForce(~,A,b,x0,jaco)

tol_rel = 0.0;
tol_abs = 0.0;
max_iter = 50;
[lambda, err, iter, flag, convergence, msg]  = fischer_newton(A, b, x0, ...
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
function [vgrad] = f_smooth(~,vgrad,N,nv)
    
vgrad = reshape(vgrad,2*N,nv);
M = ceil(N/4);
gw = gausswin(N,10);
gw = gw/sum(gw);
for nvi = 1:nv
    fc_x = vgrad(1:N,nvi);
    fc_tmp = ifft(fft(fc_x).*fft(gw));
    %    Idx = ifftshift(-N/2:N/2+1);
    %    Idx = abs(Idx) >= N/4;
    %    coeffx = fft(fc_x).*fft(gw);
    %    norm(coeffx(Idx));
    fc_x = [fc_tmp(N/2:N);fc_tmp(1:N/2-1)];
    %    coeffx = fft(fc_x);
    %    norm(coeffx(Idx));

    %coef = fft(fc_x);
    %coef(M+1:N/2+1) = 0;
    %coef(N/2+2:N-M+1) = 0;
    %fc_x = ifft(coef);

    fc_y = vgrad(N+1:2*N,nvi);
    fc_tmp = ifft(fft(fc_y).*fft(gw));
    %    coeffx = fft(fc_y).*fft(gw);
    %    norm(coeffx(Idx));
    fc_y = [fc_tmp(N/2:N);fc_tmp(1:N/2-1)];
    %coef = fft(fc_y);
    %coef(M+1:N/2+1) = 0;
    %coef(N/2+2:N-M+1) = 0;
    %fc_y = ifft(coef);

    vgrad(1:N,nvi) = fc_x;
    vgrad(N+1:2*N,nvi) = fc_y;
end
end % f_smooth

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
            vInf = [y;zeros(Np,np)];

        case 'extensional'
            vInf = [x;-y];

        case 'pipe'
            vInf = [1 - y.^2; zeros(Np,np)];
                         
        case 'taylor-green'
            vInf = [cos(x).*sin(y); -sin(x).*cos(y)];            
            
        case 'irrotational-vortex'
            vInf = [-y; x];
            
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
       M(:,k) = o.timeMatVec(I(:,k),geom,walls,true,o.explicit);
    end    
end % build_matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = build_particle_preconditioner(o)
    
    np = o.num_particles;
    Np = o.points_per_particle;
    
    P = zeros(2*Np*np + 3*np,2*Np*np + 3*np);
    
    for k = 1:np
       A = o.precop.L(:,:,k)*o.precop.U(:,:,k);
       P((k-1)*2*Np + 1: k*2*Np, (k-1)*2*Np + 1: k*2*Np) = A(1:2*Np,1:2*Np);
       P(2*Np*np + 3*(k-1)+1:2*Np*np + 3*k,(k-1)*2*Np + 1:k*2*Np) = A(end-2:end,1:end-3);
       P((k-1)*2*Np + 1:k*2*Np, 2*Np*np + 3*(k-1)+1:2*Np*np + 3*k) = A(1:end-3,end-2:end);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = build_preconditioner_matrix(o, geom, walls)

    np = o.num_particles;
    nw = o.num_walls;
    Np = o.points_per_particle;
    Nw = o.points_per_wall;
    
    Ntotal = 2*Np*np + 3*np;
    if nw > 0
        Ntotal = Ntotal + 2*Nw*nw +3*(nw-1);
    end
    
    M = zeros(Ntotal,Ntotal);
    
    % fibre-fibre
    for k = 1:np
        start_row = 1+(k-1)*2*Np;
        M(start_row:start_row+2*Np-1,start_row:start_row+2*Np-1) ...
                            =  -1/2*eye(2*Np)+o.Dp(:,:,k);
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
            [[geom.sa(:,k)'*2*pi/Np, zeros(1,Np)];...
            [zeros(1,Np), geom.sa(:,k)'*2*pi/Np]]/(2*pi);
        
        start_row = 2*Np*np+2*Nw*nw+2*np+k;
        M(start_row,start_col:start_col+2*Np-1) = ...
            [(geom.X(end/2+1:end,k)'-geom.center(2,k)).*geom.sa(:,k)'*2*pi/Np ...
            -(geom.X(1:end/2,k)'-geom.center(1,k)).*geom.sa(:,k)'*2*pi/Np]/(2*pi);
        
        start_row = 2*Np*(k-1)+1;
        start_col = 2*Np*np+2*Nw*nw+2*(k-1)+1;
        
        M(start_row:start_row+2*Np-1,start_col:start_col+1) = ...
            [[-ones(Np,1);zeros(Np,1)] ...
            [zeros(Np,1);-ones(Np,1)]];
        
        start_col = 2*Np*np + 2*Nw*nw+2*np+k;
        M(start_row:start_row+2*Np-1,start_col) = ...
            [-(geom.X(end/2+1:end,k)-geom.center(2,k)); ...
            geom.X(1:end/2,k)-geom.center(1,k)];
    end
    
    if o.confined
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
    end
end % build_preconditioner_matrix

end % methods

end % classdef


