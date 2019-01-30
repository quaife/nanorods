classdef tstep < handle
% This class defines the functions required to advance the geometry
% forward in time.

properties

tstep_order             % time stepping order
dt                      % time step size
max_dt
min_dt
num_particles           % number of particles
num_walls               % number of walls
points_per_particle     % points per particle
points_per_wall         % points per wall
Dp                      % Stokes double-layer potential for fiber-fiber interaction
Dp_old                  % Stokes double-layer potential for fiber-fiber interaction from previous time step
Dw                      % Stokes double-layer potential for wall-wall interaction
N0w                     % N0 matrix to remove rank 1 nullspace
fmm                     % flag for using the FMM
near_singular           % flag for using near-singular integration
gmres_tol               % GMRES tolerance
gmres_max_it            % maximum GMRES iterations
near_structff           % near-singular integration structure (fibre-fibre)
near_structfw           % near-singular integration structure (fibre-wall)
near_structwf           % near-singular integration structure (wall-fibre)
near_structww           % near-singular integration structure (wall-wall)
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
sedementation
torqueP0

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(options, prams, om, geom, walls, tau0, torqueP0)
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
%o.dt = 0.003125;
o.max_dt = o.dt*2;
o.min_dt = o.dt/1024;
%o.dt = o.dt/2;
o.gmres_tol = options.gmres_tol;             
o.gmres_max_it = min(300, 2*(prams.np*prams.Np + prams.nw*prams.Nw));

o.far_field = @(X,t) o.bgFlow(X, t, options); 
o.om = om;
o.tau0 = tau0;

o.num_particles = prams.np;
o.num_walls = prams.nw;
o.points_per_particle = prams.Np;
o.points_per_wall = prams.Nw;
o.debug = options.debug;
o.matvecs = 0;
o.sedementation = options.sedementation;
o.torqueP0 = torqueP0;

np = o.num_particles;
Np = o.points_per_particle;

if ~isempty(walls)
    nw = walls.n;
    Nw = walls.N;
else
    nw = 0;
    Nw = 0;
end

% CREATE CLASSES TO EVALUATE POTENTIALS ON FIBRES AND WALLS
o.potp = poten(Np, om);

if options.confined
    o.potw = poten(Nw, om);
    o.Dw = o.potw.stokesDLmatrix(walls);
    o.N0w = o.potw.stokesN0matrix(walls);
    [o.near_structww, ~] =  walls.getZone(walls,1);
else
    o.potw = [];
    o.Dw = [];
    o.N0w = [];
end

% CREATE MATRICES FOR FIBRE-FIBRE SELF INTERATIONS AND 
% WALL-WALL SELF INTERACTIONS
if ~isempty(geom.X)
    o.Dp = o.potp.stokesDLmatrix(geom);
    o.Dp_old = zeros(size(o.Dp));
else
    o.Dp = [];
    o.Dp_old = [];
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
            -(geom.X(1:end/2,k)'-geom.center(1,k)).*geom.sa(:,k)'*2*pi/Np...
            0 0 0]/(2*pi)]);
    end
    
    if o.confined
        
        % WALL-WALL PRECONDITIONER
        o.precow.L = zeros(2*Nw + 3,2*Nw + 3,nw);
        o.precow.U = zeros(2*Nw + 3,2*Nw + 3,nw);
        
        oc = curve;
        sa = walls.sa;
        [x,y] = oc.getXY(walls.X);
        [cx,cy] = oc.getXY(walls.center);
        
        for k = 1:nw
            
            if k == 1 % first wall does not need Rotlets and Stokeslets
                [o.precow.L(1:2*Nw,1:2*Nw,k),o.precow.U(1:2*Nw,1:2*Nw,k)]...
                    =lu(-1/2*eye(2*Nw)+o.Dw(:,:,k)+o.N0w(:,:,k));
                
            else
                r = [x(:,k) - cx(k), y(:,k) - cy(k)];
                rho2 = (x(:,k) - cx(k)).^2 + (y(:,k) - cy(k)).^2;
                
                col_stokes1 = [-0.5*log(rho2) + r(:,1).*r(:,1)./rho2; ...
                                r(:,2).*r(:,1)./rho2]/(4*pi);
                col_stokes2 = [r(:,2).*r(:,1)./rho2; ...
                                -0.5*log(rho2) + r(:,2).*r(:,2)./rho2]/(4*pi);
                col_rot = [r(:,2)./rho2; -r(:,1)./rho2]/(4*pi);
                
                int_stokes = sa(:,k)'/(2*Nw);
                int_rot = [(r(:,2).*sa(:,k))', -(r(:,1).*sa(:,k))']/(2*Nw);
                
                [o.precow.L(:,:,k),o.precow.U(:,:,k)] =...
                    lu([-1/2*eye(2*Nw)+o.Dw(:,:,k), col_stokes1, ...
                    col_stokes2, col_rot;...
                    int_stokes, zeros(1,Nw), -1, 0, 0;...
                    zeros(1,Nw), int_stokes, 0, -1, 0;...
                    int_rot, 0, 0, -1]);
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
                torqueP,iter,iflag,res,geomProv] = timeStep(o,geom,geomOld,walls,...
                xc,tau,uP_m1,omega_m1,first_step,forceP_old, torqueP_old,...
                etaP_old,etaW_old, time)
% timeStep(geom, walls, xc, tau) takes a single time step and returns the
% angles and velocities of the particles at the next time step as well as
% the density functions on the particles and walls, the translational and
% rotational velocities of the particles, and the net force and torque on
% all particles and walls

np = geom.n;
Np = geom.N;

if ~isempty(walls)
    nw = walls.n;
    Nw = walls.N;
else
    nw = 0;
    Nw = 0;
end

o.points_per_particle = Np;
o.num_particles = np;
o.points_per_wall = Nw;
o.num_walls = nw;

% ASSEMBLE RHS WITH NO FORCE AND TORQUE ON PARTICLES
if ~isempty(o.torqueP0)
    torqueP = o.torqueP0;
else
    torqueP = zeros(1,np);
end


forceP = zeros(2,np);
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

% CREATE MATRICES FOR FIBRE-FIBRE SELF INTERATIONS AND
% WALL-WALL SELF INTERACTIONS
if ~isempty(geom.X)
    o.Dp_old = o.Dp;
    o.Dp = o.potp.stokesDLmatrix(geom);
    %o.Dp_old = zeros(size(o.Dp));
else
    o.Dp = [];
end

% ROTATE FIBRE-FIBRE DLP AND FIBRE-FIBRE PRECONDITIONER
% dtau = tau - o.tau0;
% o.tau0 = tau;
%

o.precop.L = zeros(2*Np+3,2*Np+3,np);
o.precop.U = zeros(2*Np+3,2*Np+3,np);

%     R = spdiags([sin(dtau(i))*ones(2*Np,1), cos(dtau(i))*ones(2*Np,1)...
%                 -sin(dtau(i))*ones(2*Np,1)], [-Np, 0, Np], zeros(2*Np, 2*Np));
%
%     o.Dp(:,:,i) = R*o.Dp_old(:,:,i)*R';
%     o.Dup(:,:,i) = o.Dp(:,:,i);
%
if o.use_precond
    
    for i = 1:np
        [o.precop.L(:,:,i),o.precop.U(:,:,i)] =...
            lu([-1/2*eye(2*Np)+o.Dp(:,:,i) ...
            [-ones(Np,1);zeros(Np,1)] ...
            [zeros(Np,1);-ones(Np,1)] ...
            [-(geom.X(end/2+1:end,i)-geom.center(2,i));...
            geom.X(1:end/2,i)-geom.center(1,i)];...
            [geom.sa(:,i)'*2*pi/Np/(2*pi) zeros(1,Np) 0 0 0];
            [zeros(1,Np) geom.sa(:,i)'*2*pi/Np/(2*pi) 0 0 0];
            [(geom.X(end/2+1:end,i)'-geom.center(2,i)).*geom.sa(:,i)'*2*pi/Np ...
            -(geom.X(1:end/2,i)'-geom.center(1,i)).*geom.sa(:,i)'*2*pi/Np ...
            0 0 0]/(2*pi)]);
        
    end
    
    if o.confined
        
        % WALL-WALL PRECONDITIONER
        o.precow.L = zeros(2*Nw + 3,2*Nw + 3,nw);
        o.precow.U = zeros(2*Nw + 3,2*Nw + 3,nw);
        
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
                col_rot = [r(:,2)./rho2; -r(:,1)./rho2]/(4*pi);
                
                int_stokes = sa(:,k)'/(2*Nw);
                int_rot = [(r(:,2).*sa(:,k))', -(r(:,1).*sa(:,k))']/(2*Nw);
                
                [o.precow.L(:,:,k),o.precow.U(:,:,k)] =...
                    lu([-1/2*eye(2*Nw)+o.Dw(:,:,k), col_stokes1, ...
                    col_stokes2, col_rot;...
                    int_stokes, zeros(1,Nw), -1, 0, 0;...
                    zeros(1,Nw), int_stokes, 0, -1, 0;...
                    int_rot, 0, 0, -1]);
            end
        end
    end    
end

if explicit_time_step
    rhs = o.assembleRHS_explicit(geom, geomOld, walls, forceP, torqueP, ...
        etaP_old, 1:np, false, time);
else
    rhs = o.assembleRHS(geom, walls, forceP, torqueP, forceW, torqueW, ...
                1:np, false, time);
end

%o.gmres_max_it = 10;
% SOLVE SYSTEM USING GMRES TO GET CANDIDATE TIME STEP
if o.use_precond
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom,walls,1,1:np,1:nw,false),...
      rhs,[],o.gmres_tol,o.gmres_max_it,@o.preconditionerBD,[]);
else
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom,walls,1,1:np,1:nw,false),...
        rhs,[],o.gmres_tol, o.gmres_max_it, [], []);
end

iter = I(2);

% UNPACK SOLUTION
[etaP, etaW, uP, omega, forceW, torqueW] = o.unpackSolution(Xn,geom,...
            walls,false,true);
omega = -omega;

if np > 0
    % CREATE CANDIDATE CONFIGURATION
    if first_step || o.tstep_order == 1 % FORWARD EULER
        xc_new = xc + o.dt*uP;
        tau_new = tau + o.dt*omega;
    else %ADAMS_BASHFORTH
        xc_new = xc + (3/2)*o.dt*uP - (1/2)*o.dt*uP_m1;
        tau_new = tau + (3/2)*o.dt*omega - (1/2)*o.dt*omega_m1;
    end
    
    geomProv = copy(geom);
    geomProv.rotate(xc, tau_new - tau);
    geomProv.translate(xc_new - xc);
    
    %geomProv = capsules(o.prams, xc_new, tau_new);
    
    % RESOLVE COLLISIONS
    if o.resolve_collisions
        
        % REASSEMBLE OLD CONFIGURATION
        geomOld = copy(geom);
        
        solution.xc_old = xc;
        solution.tau_old = tau;
        solution.xc_new = xc_new;
        solution.tau_new = tau_new;
        solution.etaP = etaP;
        solution.etaW = etaW;
        solution.uP = uP;
        solution.omega = omega;
        solution.uP0 = uP;
        solution.omega0 = omega;
        solution.forceW = forceW;
        solution.torqueW = torqueW;
        solution.forceP = zeros(2,np);
        solution.torqueP = torqueP;
        solution.forceP_old = forceP_old;
        solution.torqueP_old = torqueP_old;
        solution.etaP_old = etaP_old;
        solution.etaW_old = etaW_old;
        
        solution =  o.resolveCollisions(geomOld, geomProv, walls, solution,...
            uP_m1, omega_m1, first_step,time);
        
        if o.confined
            if nw > 1
                solution.forceW = solution.forceW(:,2:end);
                solution.torqueW = solution.torqueW(2:end);
            else
                solution.forceW = [];
                solution.torqueW = [];
            end
        end
        
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
        
        geomProv = copy(geom);
        geomProv.rotate(xc, tau_new - solution.tau_old);
        geomProv.translate(xc_new - solution.xc_old);
        
    else
        if o.confined
            if nw > 1
                forceW = forceW(:,2:end);
                torqueW = torqueW(2:end);
            else
                forceW = [];
                torqueW = [];
            end
        end
    end
else
    
        if o.confined
            if nw > 1
                forceW = forceW(:,2:end);
                torqueW = torqueW(2:end);
            else
                forceW = [];
                torqueW = [];
            end
        end
        
        xc_new = xc;
        tau_new = tau;
        etaP = [];
        etaW = etaW;
        uP = [];
        omega = [];
end

end % timeStep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tx = timeMatVec(o, Xn, geom, walls, bounding_wall, ...
        particle_numbers, wall_numbers, preprocess)
% Tx = timeMatVec(Xn,geom) does a matvec for GMRES 
% if preprocess is true, solid walls are treated as particles and an
% angular and translational velocity are solved for

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

if np ~= o.num_particles || nw ~= o.num_walls
    Dp = o.Dp(:,:,particle_numbers);
    near_structff = geom.getZone([],1);
    if o.confined && nw > 0
        Dw = o.Dw(:,:,wall_numbers);
        [~,near_structfw] = geom.getZone(walls,2);
        [~,near_structwf] = walls.getZone(geom,2);
        N0w = pot_walls.stokesN0matrix(walls);
    end
else
    Dp = o.Dp;
    Dw = o.Dw;
    near_structff = o.near_structff;
    near_structfw = o.near_structfw;
    near_structwf = o.near_structwf;
    N0w = o.N0w;
end

% if nw ~= o.num_walls && nw > 0
%     Dw = o.Dw(:,:,wall_numbers);
%     %Duw = o.Duw(:,:,wall_numbers);
%     [~,near_structfw] = geom.getZone(walls,2);
%     [~,near_structwf] = walls.getZone(geom,2);
%     
%     N0w = pot_walls.stokesN0matrix(walls);
% else
%     Dw = o.Dw;
%     %Duw = o.Duw;
%     near_structfw = o.near_structfw;
%     near_structwf = o.near_structwf;
%     
%     N0w = o.N0w;
% end
    
% PREALLOCATE OUTPUT VECTORS
% velParticles : velocity of particles
% velWalls     : velocity of particles
% forceP       : net force on particles, strength of Stokeslets
% torqueP      : net torque on particles, strengh of rotlets
velParticles = zeros(2*Np,np);

% UNPACK Xn
% etaP    : denisty function on particles
% etaW    : density function on walls
% uT      : translational velocity of particles
% omega   : angular velocity of particles
% forceW  : net force on walls, strength of Stokeslets
% torqueW : net torque on walls, strength of rotlets
%if ~preprocess
    [etaP, etaW, uP, omegaP, forceW, torqueW] = unpackSolution(o, Xn,...
                geom, walls, preprocess, bounding_wall ~= 0);
% else
%     [etaP, etaW, uP, omegaP, uW, omegaW] = unpackSolution(o, Xn, geom, walls,...
%                         preprocess, bounding_wall ~= 0);
% end

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

% ADD SELF CONTRIBUTION
velParticles = velParticles + ...
    pot_particles.exactStokesDLdiag(geom, Dp, etaP);

% START OF SOURCE == PARTICLES
% START OF TARGET == PARTICLES

if ~o.explicit && np > 0
    if o.fmm
        kernel = @pot_particles.exactStokesDLfmm;
    else
        kernel = @pot_particles.exactStokesDL;
    end
    
    kernelDirect = @pot_particles.exactStokesDL;
    
    if o.near_singular        

        DLP = @(X) pot_particles.exactStokesDLdiag(geom,Dp,X) - 1/2*X;
        pp_dlp = pot_particles.nearSingInt(geom, etaP, DLP, Dp, ...
            near_structff, kernel, kernelDirect, geom, true, false);
    else
        pp_dlp = kernel(geom, etaP, Dp);
    end
else
    pp_dlp = zeros(2*Np,np);
end
% END OF TARGET == PARTICLES

% START OF TARGET == WALLS 


% IF IN PREPROCESS AND CONSIDERING COLLISION WITH BOUNDING WALL WE IGNORE
% DENSITY FUNCTION ON WALL

if o.confined && nw > 0  	
	velWalls = zeros(2*Nw,nw);
	velWalls = velWalls - 1/2*etaW;
	velWalls = velWalls + pot_walls.exactStokesDLdiag(walls, Dw, etaW);
end

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
		   near_structfw, kernel, kernelDirect, walls, false, false);
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
	   DLP = @(X) pot_particles.exactStokesDLdiag(walls, Dw, X) - 1/2*X;
	   pw_dlp = pot_particles.nearSingInt(walls, etaW, DLP, Dw, ...
		   near_structwf, kernel, kernelDirect, geom, false, false);
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
    
    kernelDirect = @pot_walls.exactStokesDL;
    
    if o.near_singular
        DLP = @(X) pot_walls.exactStokesDLdiag(walls, Dw, X) - 1/2*X;
        ww_dlp = pot_particles.nearSingInt(walls, etaW, DLP, Dw, ...
            o.near_structww, kernel, kernelDirect, walls, true, false);
    else
        ww_dlp = kernel(walls, etaW, Dw);
    end
    
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

	forceW = f_w - forceW; 
	torqueW = t_w - torqueW;
else
	if nw > 0
		[forceW, torqueW] = o.computeNetForceTorque(etaW, walls);
		p_sr = 0;
		w_sr = 0;
	else
		p_sr = 0;
		w_sr = 0;
		forceW = [];
		torqueW = [];
	end
end

% EVALUATE TOTAL VELOCITY ON PARTICLES

% ADD CONTRIBUTIONS FROM OTHER BODIES
velParticles = velParticles + p_sr + pp_dlp + pw_dlp;

% SUBTRACT VELOCITY ON SURFACE
for k = 1:np
  velParticles(1:Np,k) = velParticles(1:Np,k) - uP(1,k);
  velParticles(Np+1:end,k) = velParticles(Np+1:end,k) - uP(2,k);

  velParticles(1:Np,k) = velParticles(1:Np,k) ...
      - (geom.X(Np+1:end,k) - geom.center(2,k))*omegaP(k);
  velParticles(Np+1:end,k) = velParticles(Np+1:end,k)...
      + (geom.X(1:Np,k) - geom.center(1,k))*omegaP(k);
end

% EVALUATE VELOCITY ON WALLS
if o.confined && nw > 0
    
    velWalls = velWalls + ww_dlp + wp_dlp + w_sr;
    
    if bounding_wall ~= 0
        velWalls(:,1) = velWalls(:,1)...
            + pot_walls.exactStokesN0diag(walls, o.N0w, etaW(:,1));
    end
end


% EVALUTATE FORCES ON PARTICLES
[forceP, torqueP] = o.computeNetForceTorque(etaP, geom);

% CONSTRUCT OUTPUT VECTOR
if nw > 0    
    if bounding_wall == 0
        Tx = [velParticles(:); velWalls(:); forceP(:); torqueP(:); ...
                        forceW(:); torqueW(:)];
    else
        if nw > 1
            forceW = forceW(:,2:end);
            torqueW = torqueW(:,2:end);
        else
            forceW = [];
            torqueW = [];
        end
        
        Tx = [velParticles(:); velWalls(:); forceP(:); torqueP(:); ...
                        forceW(:); torqueW(:)];
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
function rhs = assembleRHS(o, geom, walls, forceP, torqueP, forceW, ...
                        torqueW, bodies, preprocess, time)
% Assemble right hand side
np = geom.n;
Np = geom.N;

% %%%%%% TEST %%%%%%
% torqueP = -torqueP;
% %%%%%%%%%%%%%%%%%%

if o.sedementation && ~preprocess
    forceP(2,:) = forceP(2,:) - 10;
end

if ~isempty(walls)
    nw = walls.n;
    Nw = walls.N;
else
    nw = 0;
    Nw = 0;
end

if ~preprocess
    if o.confined
        ff = o.far_field(walls.X, time);
    else
        ff = -o.far_field(geom.X, time); 
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

%if o.resolve_collisions
 
    % ADD IN CONTRIBUTIONS TO VELOCITIES FROM PARTICLE ROTLETS AND STOKESLETS
    for k = bodies
	
        % CONTRIBUTION TO PARTICLE VELOCITY
        v = o.computeRotletStokesletVelocity(geom.X, geom.center(:,k),...
            	forceP(:,k), torqueP(k));
            
        rhs(1:2*Np*np) = rhs(1:2*Np*np) - v(:);
            
    end
    
    if o.confined && ~isempty(walls)
        
        % ADD IN CONTRIBUTIONS TO VELOCITIES FROM PARTICLE ROTLETS 
        % AND STOKESLETS ON WALLS
		for k = bodies
		
			% CONTRIBUTION TO PARTICLE VELOCITY
			v = o.computeRotletStokesletVelocity(walls.X, geom.center(:,k),...
					forceP(:,k),torqueP(k));
				
			rhs(2*Np*np + 1:2*Np*np + 2*Nw*nw) = ...
                            rhs(2*Np*np + 1:2*Np*np + 2*Nw*nw) - v(:);
				
		end
    
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
%end

end % assembleRHS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhs = assembleRHS_explicit(o, geom, geomOld, walls, forceP, ...
                    torqueP, etaP_old, bodies, preprocess, time)
% Assemble right hand side

np = geom.n;
Np = geom.N;

if ~isempty(walls) && o.confined
    nw = walls.n;
    Nw = walls.N;
else
    nw = 0;
    Nw = 0;
end

if o.confined
    ff = o.far_field(walls.X, time);
else
    ff = -o.far_field(geom.X, time);
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
if ~preprocess && np > 0
    pot_particles = o.potp;
    if o.fmm
        kernel = @pot_particles.exactStokesDLfmm;
    else
        kernel = @pot_particles.exactStokesDL;
    end
    
    kernelDirect = @pot_particles.exactStokesDL;
    
    if o.near_singular
        DLP = @(X) pot_particles.exactStokesDLdiag(geom,o.Dp,X) - 1/2*X;
        pp_dlp = pot_particles.nearSingInt(geom, etaP_old, DLP, o.Dp, ...
            o.near_structff, kernel, kernelDirect, geom, true, false);
    else
        pp_dlp = kernel(geom, etaP_old, o.Dp);
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
        
        if ~preprocess && ~isempty(geomOld)
            
            % SUBTRACT OFF VELOCITY FROM ROTLETS INDUCED BY etaP_old
            [force_tmp, torque_tmp] = o.computeNetForceTorque(etaP_old, ...
                            geomOld);

            if norm(force_tmp)/length(force_tmp) > 1e-12
                
                v = o.computeRotletStokesletVelocity(geomOld.X, ...
                    geomOld.center(:,k),force_tmp(:,k),torque_tmp(k));
                
                % intra-particle interactions have already been considered
                v(:,k) = 0;
                
                rhs(1:2*Np*np) = rhs(1:2*Np*np) - v(:);
                
                if o.confined
                    % CONTRIBUTION TO WALL VELOCITY
                    v = o.computeRotletStokesletVelocity(walls.X, ...
                        geom.center(:,k),forceP(:,k),torqueP(k));
                    rhs(2*Np*np+1:2*Np*np+2*Nw*nw) = ...
                        rhs(2*Np*np+1:2*Np*np+2*Nw*nw) - v(:);
                end
            end
        end
    end
end

end % assembleRHS_explicit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [etaP, etaW, u, omega, forceW, torqueW] = unpackSolution(...
                            o, Xn, geom, walls, preprocess, bounding_wall)
% Unpack solution vector that comes from GMRES call

Np = geom.N;
np = geom.n;

if ~isempty(walls)
    Nw = walls.N;
    nw = walls.n;
  
else
    Nw = 0;
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
    etaW = [];
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
        
        stokesletStart = 2*Np*np + 2*Nw*nw + 3*np + 2*(k-2) + 1;
        stokesletEnd = stokesletStart + 1;
        rotletStart = 2*Np*np + 2*Nw*nw + 3*np + 2*(nw-1) + (k-2) + 1;
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
% x component of velocity due to the stokeslet and rotlet

LogTerm = -0.5*log(rho2)*stokeslet(2);
rorTerm = 1./rho2.*((y-cy).*(x-cx)*stokeslet(1) + ...
    (y-cy).*(y-cy)*stokeslet(2));

RotTerm = -((x-cx)./rho2)*rotlet;
vely = (LogTerm + rorTerm)/(4*pi) + RotTerm/(4*pi);
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

  force(:,k) = force(:,k)/(4*pi);
  torque(k) = torque(k)/(4*pi);  
end

end % computeNetForceTorque

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [force,torque] = computeNetForceTorqueBalanced(~, eta, geom)
% z = letsIntegrals(stokeslet,rotlet,etaM,walls) integrates the density
% function to enforce constraints on stokeslets and rotlets

N = geom.N;
n = geom.n;

force = zeros(2,n);
torque = zeros(1,n);

for k = 1:n

  force(1,k) =  sum(eta(1:N,k))*2*pi/N;
  force(2,k) =  sum(eta(N+1:end,k))*2*pi/N;

  torque(k) = sum(((geom.X(N+1:end,k) - geom.center(2,k)).*eta(1:N,k) - ...
    (geom.X(1:N,k) - geom.center(1,k)).*eta(N+1:end,k)))*2*pi/N;

  force(:,k) = force(:,k)/(4*pi);
  torque(k) = torque(k)/(4*pi);  
end

end % computeNetForceTorqueBalanced

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solution = resolveCollisions(o, geomOld, geomProv, walls, ...
                        solution, uP_m1, omega_m1, first_step,time)

Np = o.points_per_particle;
np = o.num_particles;
Nw = o.points_per_wall;
nw = o.num_walls;


upsampleFactor = 1;
nexten = 0;
c_tol = 1e-12;
minSep = o.minimum_separation;
maxIter = 1000;

% CHECK FOR COLLISION
X1 = [interpft(geomProv.X(1:Np,:),Np*upsampleFactor); ...
            interpft(geomProv.X(Np+1:end,:),Np*upsampleFactor)];
X0 = [interpft(geomOld.X(1:Np,:),Np*upsampleFactor); ...
            interpft(geomOld.X(Np+1:end,:),Np*upsampleFactor)];
geomUp = capsules([], X0);  

if o.confined
    Xw = [interpft(walls.X(1:Nw,:),Nw*upsampleFactor); ...
            interpft(walls.X(Nw+1:end,:),Nw*upsampleFactor)];
    wallsUp = capsules([], Xw);
else
    Xw = [];
    wallsUp = [];
end
    
oc = curve;
cellSizeP = 0;
cellSizeW = 0;
if np
  [~,length] = oc.geomProp(geomUp.X);
  edgelengthParticle = length/geomUp.N;
  cellSizeP = max(cellSizeP,max(edgelengthParticle));
  if o.confined
      [~,length] = oc.geomProp(wallsUp.X);
      edgelengthWall = length/wallsUp.N;
      cellSizeW = max(cellSizeW,max(edgelengthWall));
  else
	cellSizeW = cellSizeP;
  end
end

colCount = 0;

[Ns,totalnp,Xstart,Xend,totalPts] = o.preColCheck(X0, X1, Xw);

[vgrad, iv, ids, vols] = getCollision(Ns, totalnp, Xstart, Xend, minSep,...
    maxIter, totalPts, c_tol, np, nw, Np*upsampleFactor, Nw*upsampleFactor,...
    nexten, max(cellSizeP,minSep));


forceP = zeros(2,np);
forceW = zeros(2,nw);
torqueW = zeros(1,nw);

if ~isempty(o.torqueP0)
    torqueP = o.torqueP0;
else
    torqueP = zeros(1,np);
end

if iv < 0 && o.tstep_order == 2
    % compute a new candidate configuration using Forward Euler
    o.om.writeMessage('Creating candidate configuration using Forward Euler');
    
    if o.explicit
        rhs = o.assembleRHS_explicit(geomOld, walls, forceP, torqueP, ...
            solution.etaP_old, 1:np, false, time);
    else
        
        rhs = o.assembleRHS(geomOld, walls, forceP, torqueP, forceW(:,2:end), ...
            torqueW(2:end), 1:np, false, time);
    end
    
    % SOLVE SYSTEM
    if o.use_precond
        Xn = gmres(@(X) o.timeMatVec(X,geomOld,walls,1,1:np,1:nw,false),...
            rhs,[],o.gmres_tol, o.gmres_max_it,@o.preconditionerBD,[]);
    else
        Xn = gmres(@(X) o.timeMatVec(X,geomOld,walls,1,1:np,1:nw,false),...
            rhs,[],o.gmres_tol, o.gmres_max_it);
    end
    
    % UNPACK SOLUTION
    [etaP, etaW, uP, omega, forceW_new, torqueW_new] = ...
        o.unpackSolution(Xn, geomOld, walls, false, true);
    
    xc_new = solution.xc_old + o.dt*uP;
    tau_new = solution.tau_old + o.dt*omega;
    
    geomProv = capsules(o.prams, xc_new, tau_new);
    X1 = [interpft(geomProv.X(1:Np,:),Np*upsampleFactor); ...
        interpft(geomProv.X(Np+1:end,:),Np*upsampleFactor)];
    
    [Ns,totalnp,Xstart,Xend,totalPts] = o.preColCheck(X0, X1, Xw);

    [vgrad, iv, ids, vols] = getCollision(Ns, totalnp, Xstart, Xend, minSep,...
        maxIter, totalPts, c_tol, np, nw, Np*upsampleFactor, Nw*upsampleFactor,...
        nexten, max(cellSizeP,minSep));
    
    solution.xc_new = xc_new;
    solution.tau_new = tau_new;
    solution.uP = uP;
    solution.omega = omega;
    solution.forceW = forceW_new;
    solution.torqueW = torqueW_new;
    solution.forceP = forceP;
    solution.torqueP = torqueP;
    solution.etaP = etaP;
    solution.etaW = etaW;
end

n_volumes_old = max(ids);
o.om.writeMessage(['ivolume: ' num2str(iv)],'%s\n');
% RESOLVE COLLISION
while iv(end) < 0

    colCount = colCount + 1;
    o.om.writeMessage(['Collision resolving iteration: ', num2str(colCount)]);
    
    % SMOOTH VGRAD
    %vgrad(1:2*np*Np) = o.f_smooth(vgrad(1:2*np*Np),Np,np);
    
    n_volumes_new = max(ids);
    o.om.writeMessage(['Number of intersection volumes: ', num2str(n_volumes_new)])
    
    % set vgrad to normal
    nx = geomOld.xt(end/2 + 1:end,:);
    ny = -geomOld.xt(1:end/2,:);
    
    n = [nx;ny];
    n = n(:);
    n(vgrad == 0) = 0;
    vgrad = n;
    
    if colCount == 1
        n_volumes_init = max(ids);
        o.om.writeMessage(['Initial number of volumes:', num2str(n_volumes_init)]);
    end
    
    if (colCount > max(200, n_volumes_init) && n_volumes_new >= n_volumes_old && ...
                iv(end-1)/iv(end) < 2 && o.tstep_order == 1) ...
                || abs(iv(end)) > 0.1
            % halve time step if possible and volume not decreasing fast enough
            if o.dt > o.min_dt 
                o.om.writeMessage(['Reducing time step size to:', ...
                                num2str(max(o.dt/2, o.min_dt))]);
                o.dt = max(o.dt/2, o.min_dt);
                
                colCount = 1;
                
                forceP = zeros(2,np);
                forceW = zeros(2,nw);
                torqueW = zeros(1,nw);
                
                if ~isempty(o.torqueP0)
                    torqueP = o.torqueP0;
                else
                    torqueP = zeros(1,np);
                end
                
                % UPDATE PARTICLE POSTIONS AND ANGLES
                xc_new = solution.xc_old + o.dt*solution.uP0;
                tau_new = solution.tau_old + o.dt*solution.omega0;
                
                geomProv = capsules(o.prams, xc_new, tau_new);
                
                % CHECK COLLISIONS WITH UPSAMPLED GEOMETRY
                X1 = [interpft(geomProv.X(1:Np,:),Np*upsampleFactor); ...
                        interpft(geomProv.X(Np+1:end,:),Np*upsampleFactor)];
                
                [Ns,totalnp,Xstart,Xend,totalPts] = o.preColCheck(X0, X1, Xw);
                
                [vgrad, iv(end+1), ids, vols] = getCollision(Ns, totalnp, ...
                    Xstart, Xend, minSep, maxIter, totalPts, c_tol, np, ...
                    nw, Np*upsampleFactor, Nw*upsampleFactor, nexten, ...
                    max(cellSizeP,minSep));
                
            end
    else
        
        if colCount > max(100, n_volumes_init/2)  && n_volumes_new >= n_volumes_old && iv(end-1)/iv(end) < 5
            
            o.om.writeMessage('Setting vgrad to normals');
            vgrad(1:2*Np*np*upsampleFactor) = ...
                o.adjustNormal(vgrad(1:2*Np*np*upsampleFactor),geomUp,...
                edgelengthParticle);
            
            if o.confined
                vgrad(2*Np*np*upsampleFactor+1:end) = ...
                    o.adjustNormal(vgrad(2*Np*np*upsampleFactor+1:end),wallsUp,...
                    edgelengthWall);
            end
        end
    end
    
    % LINEARIZE NCP
    [A,~,ivs,~,jacoSmooth, vtoiv] = o.preprocessRigid(vgrad,ids,vols,...
        geomOld,geomUp,walls,wallsUp,torqueP);
    
    % CALCULATE FORCE AND TORQUES ON PARTICLES
    lambda = -A\(ivs/o.dt);
    fc = jacoSmooth'*lambda;
    
    % if at least one of the lambda values is negative use Fisher-Newton to
    % solve for best positive lambdas
    balance = false;
    
    if max(lambda) < 0
        o.om.writeMessage('WARNING: ALL LAMBDAS NEGATIVE, LCP CONSTRAINTS VIOLATED');
    else
        if min(lambda) < 0
            [fc, ~] = o.getColForce(A,ivs/o.dt,lambda,jacoSmooth);
            o.om.writeMessage('WARNING: NEGATIVE LAMBDA DETECTED');
            balance = true;
            %balance = false;
        end
    end
    
    fc_tmp_particles = reshape(fc(1:2*geomUp.N*np), 2*geomUp.N, np);
    
    [forceP_tmp, torqueP_tmp] = o.computeNetForceTorque(fc_tmp_particles,geomUp);

    if o.confined
        forceW_tmp = zeros(2,nw);
        torqueW_tmp = zeros(1,nw);
        bounding_wall_number = np + 1;
    else
        forceW_tmp = [];
        torqueW_tmp = [];
        bounding_wall_number = 0;
    end
    
    if balance
        [forceBalanced, torqueBalanced] = o.balance_force(geomOld, ...
            [forceP_tmp,forceW_tmp], [torqueP_tmp,torqueW_tmp], vtoiv, ...
            bounding_wall_number);
    else
        forceBalanced = [forceP_tmp,forceW_tmp];
        torqueBalanced =  [torqueP_tmp,torqueW_tmp];
    end
    
    forceP = forceP + forceBalanced(:,1:np);
    torqueP = torqueP  + torqueBalanced(1:np);
    
    % ASSEMBLE NEW RHS WITH CONTACT FORCES
    if o.explicit
        rhs = o.assembleRHS_explicit(geomOld, geomOld, walls, forceP, torqueP, ...
            solution.etaP_old, 1:np, false, time);
    else        
        rhs = o.assembleRHS(geomOld, walls, forceP, torqueP, forceW(:,2:end), ...
            torqueW(2:end), 1:np, false, time);
    end
    
    % SOLVE SYSTEM
    if o.use_precond
        Xn = gmres(@(X) o.timeMatVec(X,geomOld,walls,1,1:np,1:nw,false),...
            rhs,[],o.gmres_tol, o.gmres_max_it,@o.preconditionerBD,[]);
    else
        Xn = gmres(@(X) o.timeMatVec(X,geomOld,walls,1,1:np,1:nw,false),...
            rhs,[],o.gmres_tol, o.gmres_max_it);
    end
    
    % UNPACK SOLUTION
    [etaP, etaW, uP, omega, forceW_new, torqueW_new] = ...
            o.unpackSolution(Xn, geomOld, walls, false, true);
    
    omega = -omega;
    % UPDATE PARTICLE POSTIONS AND ANGLES
    xc_new = solution.xc_old + o.dt*uP;
    tau_new = solution.tau_old + o.dt*omega; 
    geomProv = capsules(o.prams, xc_new, tau_new);
    
    geomUp = geomProv;
    
    if o.debug
        close all
        hold on
        
        u_max_old = norm(max(abs(solution.uP), [], 2));
        f_max = norm(max(abs(forceP_tmp(:,1:o.num_particles)), [], 2));
        
        for k = 1:np
            %fill([geomUp.X(1:end/2,k);geomUp.X(1,k)],[ geomUp.X(end/2+1:end,k); geomUp.X(end/2+1,k)], 'k')
            plot([geomUp.X(1:end/2,k); geomUp.X(1,k)],[geomUp.X(end/2+1:end,k); geomUp.X(end/2+1,k)],'-ob')
            quiver(geomUp.X(1:end/2,k),geomUp.X(end/2+1:end,k), ...
                fc((k-1)*2*geomUp.N + 1:(k-1)*2*geomUp.N+geomUp.N),...
                fc((k-1)*2*geomUp.N+geomUp.N+1:(k-1)*2*geomUp.N+2*geomUp.N),'r');
            quiver(geomOld.center(1,k),geomOld.center(2,k),forceP_tmp(1,k)...
                /f_max,forceP_tmp(2,k)/f_max,1,'m');
        end
        
        
        if o.confined
            f_max = norm(max(abs(forceW_new(:,:)), [], 2));
            plot(wallsUp.X(1:end/2,:), wallsUp.X(end/2+1:end,:), '-b');
            
            for k = 2:1:nw
                startx = 2*geomUp.N*np + (k-1)*2*wallsUp.N + 1;
                starty = 2*geomUp.N*np + (k-1)*2*wallsUp.N + wallsUp.N + 1;
                
                if k > 1
                    quiver(wallsUp.center(1,k),wallsUp.center(2,k),...
                        forceW_new(1,k-1)/f_max,forceW_new(2,k-1)/f_max,1,'k');
                end
                quiver(wallsUp.X(1:end/2,k),wallsUp.X(end/2+1:end,k), ...
                    vgrad(startx:startx+wallsUp.N - 1),...
                    vgrad(starty:starty+wallsUp.N - 1),'r');
            end
        end
        
        axis equal

        pause
    end
    
    
    geomProv = capsules(o.prams, xc_new, tau_new);
    solution.xc_new = xc_new;
    solution.tau_new = tau_new;
    solution.uP = uP;
    solution.omega = omega;
    solution.forceW = forceW_new;
    solution.torqueW = torqueW_new;
    solution.forceP = forceP;
    solution.torqueP = torqueP;
    solution.etaP = etaP;
    solution.etaW = etaW;

    % CHECK COLLISIONS WITH UPSAMPLED GEOMETRY
    X1 = [interpft(geomProv.X(1:Np,:),Np*upsampleFactor); ...
        interpft(geomProv.X(Np+1:end,:),Np*upsampleFactor)];
    
    [Ns,totalnp,Xstart,Xend,totalPts] = o.preColCheck(X0, X1, Xw);
    
    [vgrad, iv(end+1), ids, vols] = getCollision(Ns, totalnp, Xstart, ...
        Xend, minSep, maxIter, totalPts, c_tol, np, nw, ...
        Np*upsampleFactor, Nw*upsampleFactor,nexten, max(cellSizeP,minSep));
    
    n_volumes_old = n_volumes_new;
    
    o.om.writeMessage(['ivolume: ' num2str(iv(end))],'%s\n');
end

% INCREASE TIME STEP SIZE IF DESIRED
% if colCount <= 3 && o.dt < o.max_dt && o.tstep_order == 1
%     o.om.writeMessage(['Increasing time step size to:', ...
%         num2str(min(o.max_dt, 1.2*o.dt))]);
%     o.dt = min(o.max_dt, 1.5*o.dt);
% end


end % resolveCollision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,jaco,ivs,listnv,jacoSmooth, vtoiv] = preprocessRigid(o,vgrad,ids,...
        vols, geomOld, geomUp, walls, wallsUp,torqueP)

Np = geomUp.N;
np = geomUp.n;

if ~isempty(wallsUp)
nw = wallsUp.n;
Nw = wallsUp.N;
else
nw = 0;
Nw = 0;
end

nivs = max(ids);
A = zeros(nivs,nivs);
ivs = zeros(nivs,1);
%vtoiv = zeros(np + nw,nivs);
vtoiv = zeros(np, nivs);

jI = [];
jJ = [];
jV = [];
jVr = [];
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
jacoSmooth = jaco;

% COMPUTE AND BALANCE FORCES

f = jacoSmooth'*ones(nivs,1);
fp = reshape(f(1:2*Np*np), 2*Np, np);
[forceP_orig,torqueP_orig] = o.computeNetForceTorque(fp,geomOld);
%[forceP_all,torqueP_all] = o.computeNetForceTorqueBalanced(fp,geomUp);

if o.confined
%   fw = reshape(f(2*Np*np+1:end), 2*Nw, nw);
%   [forceW_orig,torqueW_orig] = o.computeNetForceTorque(fw,wallsUp);
    forceW_orig = zeros(2,nw);
    torqueW_orig = zeros(1,nw);
    bounding_wall_number = np + 1;
else
    forceW_orig = [];
    torqueW_orig = [];
    bounding_wall_number = 0;
end


% [forceBalanced, torqueBalanced] = o.balance_force(geomOld, ...
%             [forceP_orig,forceW_orig], [torqueP_orig, torqueW_orig],...
%             vtoiv, bounding_wall_number);

forceBalanced = [forceP_orig,forceW_orig];
torqueBalanced =  [torqueP_orig, torqueW_orig];

forceP_all = forceBalanced(:,1:np);
torqueP_all = torqueBalanced(1:np);

forceW_all = zeros(2,nw);
torqueW_all = zeros(nw,1);

%C = o.determine_clusters(vtoiv, bounding_wall_number);

% LOOP OVER INTERSECTION VOLUMES
for i = 1:nivs
    
    % FIND BODIES IN THIS VOLUME
    S = find(vtoiv(:,i)~=0); 
    
    forceP = forceP_all(:,S(S<=np)');
    torqueP = torqueP_all(S(S<=np)');
    
    if o.confined
        forceW = forceW_all(:,S(S>np) - np);
        torqueW = torqueW_all(S(S>np) - np);
    else
        forceW = [];
        torqueW = [];
    end
    
    geomTmp = capsules([], geomOld.X(:,S(S<=np)'));

    if max(S) > np
        wallsTmp = capsules([], walls.X(:, S(S > np) - np));
    else
        wallsTmp = [];
    end
 
    if max(S == np+1) % collision with outer wall
        bounding_wall = 1;
        
        if wallsTmp.n == 1
            forceW = [];
            torqueW = [];
        end
    else
        bounding_wall = 0;
    end
    
    if o.confined
        for c = 1:size(C,2)
            if max(C{c}' == S(1))
                cluster = C{c}';
                n_walls = length(cluster(cluster > np));
            end
        end
        
        if n_walls > 0 %if collision involves wall, scale force on particle
            forceP = forceP/(2*n_walls);
            torqueP = torqueP/(2*n_walls);
        end
    end
    
    % ASSEMBLE NEW RHS WITH CONTACT FORCES
    if o.explicit 
        rhs = o.assembleRHS_explicit(geomTmp, geomOld, wallsTmp, forceP, torqueP, ...
                [], 1:length(S(S<=np)), true, 0);
    else
        rhs = o.assembleRHS(geomTmp, wallsTmp, forceP, torqueP, forceW,...
                        torqueW, 1:length(S(S<=np)), true, 0);  
    end
    
    % SOLVE SYSTEM WITH FAR FIELD NEGLECTED   
    Xn = gmres(@(X) o.timeMatVec(X,geomTmp,wallsTmp,bounding_wall,...
                S(S<=np),S(S > np) - np,true),rhs, [],o.gmres_tol,length(rhs));
    
    % COMPUTE VELOCITY OF EACH POINT ON ALL RIGID PARTICLES
    [~, ~, uP, omega, ~, ~] = o.unpackSolution(Xn, geomTmp, wallsTmp, ...
                true, bounding_wall == 1);
    
    for k = 1:geomTmp.n
        b = [uP(1,k)*ones(geomUp.N,1); uP(2,k)*ones(geomUp.N,1)] + ... %translational
            omega(k)*[geomUp.X(end/2+1:end,S(k)) - geomTmp.center(2,k); ...% + rotational
            -(geomUp.X(1:end/2,S(k)) - geomTmp.center(1,k))];
        
        j = S(k);
        SS = find(vtoiv(j,:)~=0);
        for l = 1:numel(SS)
            A(SS(l),i) = A(SS(l),i) + dot(jaco(SS(l),1+(j-1)...
                *2*geomUp.N:2*geomUp.N+(j-1)*2*geomUp.N),b);
        end
    end
    
    if o.debug
        close all
        hold on
        for k = 1:geomTmp.n
            plot(geomUp.X(1:end/2,S(k)),geomUp.X(end/2+1:end,S(k)),'b');
            quiver(geomUp.X(1:end/2,S(k)),geomUp.X(end/2+1:end,S(k)), ...
                vgrad((S(k)-1)*2*Np + 1:(S(k)-1)*2*Np+Np),...
                vgrad((S(k)-1)*2*Np+Np+1:(S(k)-1)*2*Np+2*Np),'r');
            quiver(geomTmp.center(1,k),geomTmp.center(2,k),...
                    uP(1,k),uP(2,k),1000,'g');
            quiver(geomTmp.center(1,k),geomTmp.center(2,k),...
                forceP(1,k),forceP(2,k),1000,'m');

            b = [uP(1,k)*ones(geomUp.N,1); uP(2,k)*ones(geomUp.N,1)] + ... %translational
            omega(k)*[geomUp.X(end/2+1:end,S(k)) - geomTmp.center(2,k); ...% + rotational
            -(geomUp.X(1:end/2,S(k)) - geomTmp.center(1,k))];
           % b = o.computeRotletStokesletVelocity(geomTmp.X(:,k),...
           % geomTmp.center(:,1),forceP(:,1),0);
            quiver(geomUp.X(1:end/2,S(k)),geomUp.X(end/2+1:end,S(k)), ...
                b(1:end/2),b(end/2+1:end),'k');
        end
        
        if o.confined && ~isempty(wallsTmp)
            for k = 1:wallsTmp.n
                Sw = S(S > np);
                k1 = Sw(k) - np;
                
                startx = 2*geomOld.N*np + (k1-1)*2*walls.N + 1;
                starty = 2*geomOld.N*np + (k1-1)*2*walls.N + walls.N + 1;
                plot(wallsTmp.X(1:end/2,k),wallsTmp.X(end/2+1:end,k),'-bo');
                quiver(wallsTmp.X(1:end/2,k),wallsTmp.X(end/2+1:end,k), ...
                    f(startx:startx + wallsTmp.N - 1),...
                    f(starty:starty + wallsTmp.N - 1),'r');
            end
        end
        
        axis equal
        pause
    end
end

end % preprocessRigid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F, L] = balance_force(o, geom, F, L, vtoiv, bounding_wall_number)
    
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

C = o.determine_clusters(vtoiv, bounding_wall_number);

for ic = 1:length(C)
    
    % identify which particle is in the most volumes (i.e. has the most
    % neighbors)
    
    C_tmp = C{ic};
    
    % don't balance forces for clusters containing a solid wall
    if ~max(C_tmp > geom.n)
        
        % remove particles that actually have zero net force and torque
        i_remove = [];
        for i = 1:length(C_tmp)
            if L(C_tmp(i)) == 0 && C_tmp(i) <= geom.n
                i_remove = [i_remove, i];
            end
        end
        
        C_tmp(i_remove) = [];
        
        if ~isempty(C_tmp(C_tmp <= geom.n))
            if size(vtoiv,2) > 1
                max_neighbors = max(sum(vtoiv(C_tmp,:)'));
                imax_neighbors = find(sum(vtoiv(C_tmp,:)') == max_neighbors);
            else
                max_neighbors = max(vtoiv);
                imax_neighbors = find(vtoiv(C_tmp,:)' == max_neighbors);
            end
            
            max_cluster = C_tmp(imax_neighbors);
            target_forces = F(:,max_cluster);
            [~, itarget] = max(sum(target_forces.^2));
            target_index = max_cluster(itarget);
            
            % total force on this particle will be negative sum of forces on all
            % other particles in cluser
            if length(C_tmp(C_tmp <= geom.n)) > 1 && max(C_tmp) <= geom.n
                total_force = [0;0];
                for k = C_tmp'
                    if k ~= target_index
                        total_force = total_force + F(:,k);
                    end
                end
                
                f_orig = F(:,target_index);
                F(:,target_index) = -total_force;
                
                L(target_index) = L(target_index)...
                                *norm(F(:,target_index))/norm(f_orig);
            end
        end
    end
end

end % balance_force

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C, volumes] = add_particle_to_cluster(o, C, S, k, volumes, ...
                    bounding_wall_number)
    
n_vols = length(S);

for i = 1:n_vols
    if ~max(i == volumes)
        if max(S{i} == k) && k ~= bounding_wall_number
            volumes(end+1) = i;
            C = [C; S{i}];

            for k1 = S{i}'
                [C, volumes] = o.add_particle_to_cluster(C, S, k1, ...
                    volumes, bounding_wall_number);
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = determine_clusters(o, vtoiv, bounding_wall_number)
    
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
             [C{end+1}, already_considered] = o.add_particle_to_cluster([], S,...
                 k, already_considered, bounding_wall_number);
        end
    end
end

C = C(~cellfun('isempty', C));

for i = 1:length(C)
   C{i} = unique(C{i}); 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vgrad_new = correct_for_rotation(~, vgrad, geom, k)
    
    A = [[ones(geom.N,1);zeros(geom.N,1)], [zeros(geom.N,1);ones(geom.N,1)],...
                [(geom.X(end/2+1:end,k) - geom.center(2,k))...
                    ;-(geom.X(1:end/2,k) - geom.center(1,k))]];
                
    U = A\vgrad';
    
    vgrad_new = A*U;
                
    vgrad_new = vgrad_new(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vgrad = adjustNormal(~,vgrad,geom,edgelength)
% use the normal of Xn or X?
N = geom.N;
n = geom.n;

vgrad = reshape(vgrad,2*N,n);

for k = 1:n
    for i = 1:N
        if(vgrad(i,k)==0 && vgrad(N+i,k)==0)
            continue;
        end
        
        % in rare case if space-time gradient is too small,
        % collision resolving may get stuck
        % use surface normal instead
        
        tangent = geom.xt;
        oc = curve;
        [tanx,tany] = oc.getXY(tangent);
        nx = tany;
        ny = -tanx;
        normal = [nx;ny];
        
        bnnorm = edgelength(k);
        vgrad(i,k) = normal(i,k)*bnnorm;
        vgrad(N+i,k) = normal(N+i,k)*bnnorm;
    end
end

vgrad = vgrad(:);

end % adjustNormal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fc,lambda] = getColForce(~,A,b,x0,jaco)

tol_rel = 0.0;
tol_abs = 0.0;
max_iter = 50;
[lambda, ~, ~, ~, ~, ~]  = fischer_newton(A, b, x0, ...
                max_iter, tol_rel, tol_abs, 'perturbation', false);
            
fc = jaco'*lambda;
end % getColForce

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ns,totalnv,Xstart,Xend,totalPts] = preColCheck(~,X0,X1,Xw)
    
[Np, np] = size(X0);
[Nw, nw] = size(Xw);

Np = Np/2;
Nw = Nw/2;

Ns = ones(np,1)*Np;
totalnv = np;

% upsample  positions
Xv0 = reshape(X0,Np,2*np);
Xv1 = reshape(X1,Np,2*np);

Xstart = Xv0(:);
Xend = Xv1(:);

if nw > 0
  Ns = [Ns;ones(nw,1)*Nw];
  totalnv = totalnv + nw;
  
  Xb = reshape(Xw,Nw,2*nw);
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
function vInf = bgFlow(~,X,time,options)
    
N = size(X,1)/2;
n = size(X,2);
oc = curve;

[x,y] = oc.getXY(X);

if options.confined
    switch options.far_field
        case 'constant'
            vInf = [ones(N,n);zeros(N,n)];

        case 'circle'
            vInf = [-y(:,1); x(:,1)];
            
        case 'couette'

            vInf = zeros(2*N,1);
            vInf = [[-options.couette_speed*y(:,1);...
                            options.couette_speed*x(:,1)], vInf];
            
        case 'pipe' 
            if time < 5
                vInf = [[1 - y(:,1).^2; zeros(N,1)], zeros(2*N,1)];
            else
                vInf = -[[1 - y(:,1).^2; zeros(N,1)], zeros(2*N,1)];
            end
            
        case 'bounded_shear'            
            if n == 1
                vInf = [-y(:,1); zeros(N,1)];
            else if n == 2
                    vInf = [[-y(:,1);zeros(2*N,1)],[-y(:,2);zeros(2*N,1)]];
                end
            end
         
        case 'bounded_extensional'            
            if n == 1
                vInf = [x(:,1); -y(:,1)];
            end
           
        case 'bounded_eigenvalues'
            vInf = zeros(2*N,3);
            vInf = [[-y(:,1);x(:,1)], vInf];
              %vInf = [-y(:,1);x(:,1)];
              
        case 'injection_plate'
            
            vInf = zeros(2*N,n);
            Q = 2;
            r_in = 0.5;
            r_out = 10;
            
            factor = 1;
            if time > 4
                factor = -1;
            else
                factor = 1;
            end
            
            % outer wall
            vInf(:,1) = factor*Q*r_in/r_out*[x(:,1); y(:,1)];
            
            % inner wall
            vInf(:,2) = factor*Q*[x(:,2); y(:,2)];
            
    
        otherwise
            vInf = [zeros(N,n);zeros(N,n)];
    end
else
    switch options.far_field
        case 'shear'
            vInf = -1*[y;zeros(N,n)];

        case 'shear_reverse'
            if time < 10
                vInf = [3*y;zeros(N,n)];
            else
                vInf = -[3*y;zeros(N,n)];
            end
            
        case 'extensional'
            vInf = [x;-y];

        case 'pipe'
            vInf = [1 - y.^2; zeros(N,n)];
                         
        case 'taylor-green'
            vInf = [cos(x).*sin(y); -sin(x).*cos(y)];            
            
        case 'irrotational-vortex'
            vInf = [-y; x];
        
	case 'none'
		vInf = [zeros(N,n);zeros(N,n)];
 
        otherwise
            vInf = [zeros(N,n);zeros(N,n)];
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
        if geom.n > 0
            p = 1:geom.n;
        else
            p = [];
        end
        
        if ~isempty(walls) > 0
            w = 1:walls.n;
        else
            w = [];
        end
        
       M(:,k) = o.timeMatVec(I(:,k),geom, walls, 1, p, w, false);
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
            col_rot = [r(:,2)./rho2; -r(:,1)./rho2]/(4*pi);
            
            int_stokes = 2*pi*sa(:,k+1)'/(2*Nw);
            int_rot = 2*pi*[(y(:,k+1).*sa(:,k+1))', -(x(:,k+1).*sa(:,k+1))']/(2*Nw);
            
            start_row = 2*Np*np+2*Nw*nw+3*np+1+3*(k-1);
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

function [M_permute, P_permute] = permute_matrix(o, M)
    M_permute = zeros(size(M));
    P_permute = zeros(size(M));
    
    np = o.num_particles;
    nw = o.num_walls;
    Np = o.points_per_particle;
    Nw = o.points_per_wall;
    
    start_stokes = 2*Np*np + 2*Nw*nw + 1;
    
    for k = 1:np
        for j = 1:np
            
            % extract blocks from M
            start_hydro = 2*Np*(k-1) + 1;
            start_hydro_target = 2*Np*(j-1) + 1;
            
            % extract hydrodynamic interactions
            block = M(start_hydro:start_hydro+2*Np-1,...
                start_hydro_target:start_hydro_target+2*Np-1);
            
            if k == j
                start_force_row = start_stokes + 2*(k-1);
                start_force_column = start_stokes + 2*(j-1);
                force_block_hor = M(start_force_row:start_force_row+1,start_hydro:start_hydro+2*Np-1);
                force_block_vert = M(start_hydro_target:start_hydro_target+2*Np-1, start_force_column:start_force_column+1);
                
                torque_row = M(start_stokes + 2*np + k-1,start_hydro:start_hydro+2*Np-1);
                torque_column = M(start_hydro_target:start_hydro_target+2*Np-1,start_stokes + 2*np + j-1);
                block = [[block, force_block_vert,torque_column]; [[force_block_hor; torque_row],zeros(3)]];
                
                start_permute_row = (k-1)*2*Np +3*(k-1) + 1;
                start_permute_column = (j-1)*2*Np +3*(j-1) + 1;
                M_permute(start_permute_row : start_permute_row+2*Np+2, start_permute_column : start_permute_column+2*Np+2) = block;
                P_permute(start_permute_row : start_permute_row+2*Np+2, start_permute_column : start_permute_column+2*Np+2) = block;
            else
                
                block = [[block, zeros(2*Np,3)]; [zeros(3,2*Np), zeros(3)]];
                
                start_permute_row = (k-1)*2*Np +3*(k-1) + 1;            
                start_permute_column = (j-1)*2*Np +3*(j-1) + 1;
                
                
                M_permute(start_permute_row : start_permute_row+2*Np+2, start_permute_column : start_permute_column+2*Np+2) = block;
            end
        end
    end
    
    start_walls = 2*Np*np+1;
    start_stokes = 2*Np*np + 2*Nw*nw + 3*np + 1;
    for k = 1:nw
        for j = 1:nw
            
            % extract blocks from M
            start_hydro = start_walls + 2*Nw*(k-1);
            start_hydro_target = start_walls + 2*Nw*(j-1);
            
            % extract hydrodynamic interactions
            block = M(start_hydro:start_hydro+2*Nw-1,...
                start_hydro_target:start_hydro_target+2*Nw-1);
            
            if k == j
                if k > 1 && j > 1
                    start_force_row = start_stokes + 2*(k-2);
                    start_force_column = start_stokes + 2*(j-2);
                    force_block_hor = M(start_force_row:start_force_row+1,start_hydro:start_hydro+2*Nw-1);
                    force_block_vert = M(start_hydro_target:start_hydro_target+2*Nw-1, start_force_column:start_force_column+1);
                    
                    torque_row = M(start_stokes + 2*(nw-1) + k-2,start_hydro:start_hydro+2*Nw-1);
                    torque_column = M(start_hydro_target:start_hydro_target+2*Nw-1,start_stokes + 2*(nw-1) + j-2);
                    block = [[block, force_block_vert,torque_column]; [[force_block_hor; torque_row],-eye(3)]];
                    
                    start_permute_row = 2*Np*np+ (k-1)*2*Nw +3*(k-2) + 1;
                    start_permute_column = 2*Np*np+ (j-1)*2*Nw +3*(j-2) + 1;
                    M_permute(start_permute_row : start_permute_row+2*Nw+2, start_permute_column : start_permute_column+2*Nw+2) = block;
                    P_permute(start_permute_row : start_permute_row+2*Nw+2, start_permute_column : start_permute_column+2*Nw+2) = block;
                else
                    
                    start_permute_row = 2*Np*np+ (k-1)*2*Nw + 1;
                    start_permute_column = 2*Np*np+ (j-1)*2*Nw + 1;
                    
                    M_permute(start_permute_row : start_permute_row+2*Nw-1, start_permute_column : start_permute_column+2*Nw-1) = block;
                    P_permute(start_permute_row : start_permute_row+2*Nw-1, start_permute_column : start_permute_column+2*Nw-1) = block;
                end
            else
                block = [[block, zeros(2*Nw,3)]; [zeros(3,2*Nw), zeros(3)]];
                if k > 1
                    start_permute_row = 2*Np*np+ (k-1)*2*Nw +3*(k-2) + 1;
                else
                    start_permute_row = 2*Np*np+ (k-1)*2*Nw +1;
                end
                
                if j > 1
                    start_permute_column = 2*Np*np+ (j-1)*2*Nw +3*(j-2) + 1;
                else
                    start_permute_column = 2*Np*np+ (j-1)*2*Nw + 1;
                end
                
                M_permute(start_permute_row : start_permute_row+2*Nw+2, start_permute_column : start_permute_column+2*Nw+2) = block;
            end
        end
    end
end
end % methods

end % classdef


