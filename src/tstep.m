classdef tstep < handle
% This class defines the functions required to advance the geometry
% forward in time.  Routines that we may need to add are different
% integrators such as semi-implicit, SDC, Runge-Kutta, adaptivity

properties

tstep_order          % time stepping order
dt             % time step size
Dp             % Stokes double-layer potential for fiber-fiber interaction
Dw             % Stokes double-layer potential for wall-wall interaction
Dup           % Upsampled Stokes double-layer potential matrix for fibers
Duw           % Upsampled Stokes double-layer potential matrix for walls
N0w            % N0 matrix to remove rank 1 nullspace
rhs            % Right hand side
fmm            % flag for using the FMM
near_singular  % flag for using near-singular integration
gmres_tol     % GMRES tolerance
near_structff % near-singular integration structure (fibre-fibre)
near_structfw % near-singular integration structure (fibre-wall)
near_structwf % near-singular integration structure (wall-fibre)
far_field     % background flow
use_precond     % use a block-diagonal preconditioner
precop       % block-diagonal preconditioner for fibres
precow       % block-diagonal preconditioner for walls
potp     % class for fiber layer potentials
potw      % class for wall layer potentials
profile      % flag to time certain parts of code
om           % monitor class
confined     % flag to indicate whether flow is bounded
tau0         % initial angles of fibres
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

o.dt = prams.T/prams.number_steps;                  
o.gmres_tol = options.gmres_tol;             

o.far_field = @(X) o.bgFlow(X,options); 
o.om = om;
o.tau0 = tau0;

Np = prams.Np;
Nw = prams.Nw;
np = prams.np;
nw = prams.nw;

% CREATE CLASSES TO EVALUATE POTENTIALS ON FIBRES AND WALLS
o.potp = poten(prams.Np, om);

if options.confined
    o.potw = poten(prams.Nw, om);
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
            [geom.X(end/2+1:end,k)-geom.center(2,k); -geom.X(1:end/2,k)+geom.center(1,k)];...
            [-geom.sa(:,k)'*2*pi/Np zeros(1,Np) 0 0 0];
            [zeros(1,Np) -geom.sa(:,k)'*2*pi/Np 0 0 0];
            [(-geom.X(end/2+1:end,k)'+geom.center(2,k)).*geom.sa(:,k)'*2*pi/Np ...
            (geom.X(1:end/2,k)'-geom.center(1,k)).*geom.sa(:,k)'*2*pi/Np 0 0 0]]);
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
                
                col_stokes1 = [-0.5*log(rho2) + r(:,1).*r(:,1)./rho2; r(:,2).*r(:,1)./rho2]/(4*pi);
                col_stokes2 = [r(:,2).*r(:,1)./rho2; -0.5*log(rho2) + r(:,2).*r(:,2)./rho2]/(4*pi);
                col_rot = [r(:,2)./rho2; -r(:,1)./rho2];
                
                int_stokes = 2*pi*sa(:,k)'/Nw;
                int_rot = 2*pi*[(y(:,k).*sa(:,k))', -(x(:,k).*sa(:,k))']/Nw;
                
                [o.precow.L(:,:,k),o.precow.U(:,:,k)] =...
                    lu([-1/2*eye(2*Nw)+o.Dw(:,:,k), col_stokes1, col_stokes2, col_rot;...
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

% CREATE RIGHT HAND SIDE
if options.confined
    rhs = o.far_field(walls.X);
    o.rhs = [zeros(2*Np*np,1); rhs(:); zeros(3*np,1);zeros(3*(nw-1),1)];
else
    rhs = o.far_field(geom.X);
    o.rhs = [-rhs(:); zeros(2*Nw*nw,1); zeros(3*np,1)];
end

end % constructor: tstep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [etaF,etaW,Up,wp,stokes,rot,iter,iflag,res] = ...
                        timeStep(o,geom,tau,walls,options,prams)
% [X,iter,iflag] = timeStep(Xstore) takes the current configuration
% in Xstore (can be a three-dimensional array corresponding to  previous
% time steps if doing multistep) and returns a new shape X, the number
% of GMRES iterations, and the GMRES flag 


if ~o.confined % update RHS
    o.far_field = @(X) o.bgFlow(X,options); 
    
    Nw = prams.Nw;
    np = prams.np;
    nw = prams.nw;
    
    rhs = o.far_field(geom.X);
    o.rhs = [-rhs(:); zeros(2*Nw*nw,1); zeros(3*np,1)];
end

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
Np = geom.N;
np = geom.n;

if options.confined
    Nw = walls.N;
    nw = walls.n;
else
    Nw = 0;
    nw = 0;
end

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
            [geom.X(end/2+1:end,i)-geom.center(2,i); -geom.X(1:end/2,i)+geom.center(1,i)];...
            [-geom.sa(:,i)'*2*pi/Np zeros(1,Np) 0 0 0];
            [zeros(1,Np) -geom.sa(:,i)'*2*pi/Np 0 0 0];
            [(-geom.X(end/2+1:end,i)'+geom.center(2,i)).*geom.sa(:,i)'*2*pi/Np ...
            (geom.X(1:end/2,i)'-geom.center(1,i)).*geom.sa(:,i)'*2*pi/Np 0 0 0]]);
    end    
end

% SOLVE SYSTEM USING GMRES
maxit = 2*(Np*np + Nw*nw);% should be a lot lower than this
%o.rhs  = [ones(2*nv*N,1);2*ones(2*Nbd*nbd,1);3*ones(2*nv,1);4*ones(nv,1);5*ones(2*(nbd-1),1);6*ones(nbd-1,1)];
% 
% clf;
% semilogy(abs(o.timeMatVec(o.preconditionerBD(o.rhs), geom, walls)-o.rhs), '.');

if o.use_precond
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom,walls),o.rhs,[],o.gmres_tol,...
      maxit,@o.preconditionerBD,[], o.rhs);
else
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom,walls),o.rhs,[],o.gmres_tol,...
      maxit);
end

iter = I(2);

% REORGANIZE COLUMN VECTOR INTO MATRIX
% EXTRACT DENSITY FUNCITONS ON FIBRES AND WALLS
% each column of etaF corresponds to the density function of a rigid body
etaF = zeros(2*Np,np);
for k = 1:np
  etaF(:,k) = Xn((k-1)*2*Np+1:k*2*Np);
end

% each column of etaW corresponds to the density function of a solid wall
if options.confined
    etaW = zeros(2*Nw, nw);
    for k = 1:nw
       etaW(:,k) = Xn(2*Np*np+1+(k-1)*2*Nw:2*Np*np+k*2*Nw);
    end
else
    etaW = 0;
end

% EXTRACT TRANSLATIONAL AND ROTATIONAL VELOCITIES
Up = zeros(2,np);
for k = 1:np
  Up(:,k) = Xn(2*Np*np+2*Nw*nw+(k-1)*2+1:2*Np*np+2*Nw*nw+k*2);
end

wp = zeros(1,np);
for k = 1:np
  wp(k) = Xn(2*Np*np+2*Nw*nw+2*np+k);
end

if options.confined
    stokes = zeros(2,nw-1);
    for k = 1:nw-1
       stokes(:,k) = Xn(2*Np*np+2*Nw*nw+3*np+1+2*(k-1): 2*Np*np+2*Nw*nw+3*np+2*k);
    end

    rot = zeros(1,nw-1);
    for k = 1:nw-1
       rot(k) = Xn(2*Np*np+2*Nw*nw+3*np+2*(nw-1)+k);
    end
else
    stokes = 0;
    rot = 0;
end

end % timeStep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tx = timeMatVec(o,Xn,geom,walls)
% Tx = timeMatVec(Xn,geom) does a matvec for GMRES 
% Xn is a state vector that contains the following in order:
% fiber1: eta_x, eta_y fiber2: eta_x, eta_y ... fibernv: eta_x, eta_y
% wall1 : xi_x, xi_y wall2: xi_x, xi_y ... wallnbd: xi_x, xi_y
% fiber1: u, v fiber2: u, v ... fibernv: u, v
% fiber1: omega fiber2: omega ... fibernv: omega
% wall2 : stokeslets ... wallnbd stokeslets 
% wall2 : rotlet ... wallnbd rotlet
if o.profile
    tMatvec = tic;
end

Np = geom.N;   % points per body
np = geom.n; % number of bodies

if o.confined
    Nw = walls.N;   % points per wall
    nw = walls.n; % number of wallls
else
    Nw = 0;
    nw = 0;
end

pot_particles = o.potp;
pot_walls = o.potw;

% Output of Tx that corresponds to the velocity of the fibers
valFibers = zeros(2*Np,np);
% Output of Tx that corresponds to the velocity of the walls
valWalls = zeros(2*Nw,nw);
% output of Tx that corresponds to force on fibres
valForce = zeros(2,np);
% output of Tx that corresponds to torque on fibres
valTorque = zeros(np,1);

% BEGIN FORMATTING UNKNOWN VECTOR

% EXTRACT DENSITY FUNCTION FOR FIBRES
etaF = zeros(2*Np,np);
for k = 1:np
  etaF(:,k) = Xn((k-1)*2*Np+1:k*2*Np);
end

% EXTRACT DENSITY FUNCTION FOR WALLS
etaW = zeros(2*Nw,nw);
if o.confined
    for k = 1:nw
      etaW(:,k) = Xn(2*Np*np+1+(k-1)*2*Nw:2*Np*np+k*2*Nw);
    end
end

% EXTRACT TRANSLATIONAL AND ROTATIONAL VELOCITIES OF FIBRES
Up = zeros(2,np);
for k = 1:np
  Up(:,k) = Xn(2*Np*np+2*Nw*nw+1+2*(k-1):2*Np*np+2*Nw*nw+2*k);
end
wp = zeros(1,np);
for k = 1:np
  wp(k) = Xn(2*Np*np+2*Nw*nw+2*np+k);
end

if o.confined
    % EXTRACT ROTLETS AND STOKESLETS
    lambda = zeros(2,nw-1);
    for k = 1:nw-1
       lambda(:,k) = Xn(2*Np*np+2*Nw*nw+3*np+1+2*(k-1):2*Np*np+2*Nw*nw+3*np+2*k); 
    end

    xi=zeros(1,nw-1);
    for k=1:nw-1
        xi(k) = Xn(2*Np*np+2*Nw*nw+3*np+2*(nw-1)+k);
    end    
else
    lambda = zeros(2,1);
    xi = 0;
end

% END FORMATTING UNKNOWN VECTOR

% ADD JUMP IN DLP
valFibers = valFibers - 1/2*etaF;
if o.confined
   valWalls = valWalls - 1/2*etaW;
end 

% ADD SELF CONTRIBUTION
valFibers = valFibers + pot_particles.exactStokesDLdiag(geom, o.Dp, etaF);
if o.confined
    valWalls = valWalls + pot_walls.exactStokesDLdiag(walls, o.Dw, etaW);
end

% START OF SOURCE == FIBRES
% START OF TARGET == FIBRES
if o.fmm
  kernel = @pot_particles.exactStokesDLfmm;
else
  kernel = @pot_particles.exactStokesDL;
end

kernelDirect = @pot_particles.exactStokesDL;

if o.near_singular 
  DLP = @(X) pot_particles.exactStokesDLdiag(geom,o.Dp,X) - 1/2*X;
  ffdlp = pot_particles.nearSingInt(geom, etaF, DLP, o.Dup, o.near_structff ,kernel,...
        kernelDirect, geom, true, false);
else
  ffdlp = kernel(geom, etaF, o.Dp);
end
% END OF TARGET == FIBRES

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
       wfdlp = pot_walls.nearSingInt(geom, etaF, DLP, o.Dup, o.near_structfw, ...
           kernel, kernelDirect, walls, false, false);
   else
       wfdlp = kernel(geom, etaF);
   end
else
    wfdlp = zeros(2*Nw,nw);
end
% END OF TARGET == WALLS

% START OF SOURCE == WALLS
% START OF TARGET == FIBERS
if o.confined && np > 0
   
   if o.fmm
        kernel = @pot_particles.exactStokesDLfmm;       
   else
       kernel = @pot_particles.exactStokesDL;
   end
   
   kernelDirect = @pot_particles.exactStokesDL;
   
   if o.near_singular
       DLP = @(X) pot_particles.exactStokesDLdiag(walls, o.Dw, X) - 1/2*X;
       fwdlp = pot_particles.nearSingInt(walls, etaW, DLP, o.Duw, o.near_structwf, ...
           kernel, kernelDirect, geom, false, false);
   else
       fwdlp = kernel(walls, etaW);
   end
else
    fwdlp = zeros(2*Np,np);
end
% END OF TARGET == FIBRES

% START OF TARGET == WALLS
if o.confined
    
    if o.fmm
        kernel = @pot_walls.exactStokesDLfmm;
    else
        kernel = @pot_walls.exactStokesDL;
    end
    
    wwdlp = kernel(walls, etaW, o.Dw);
    
else
    wwdlp = zeros(2*Nw,nw);
end
% END OF TARGET == WALLS
% END OF SOURCE == WALLS

if (o.confined && nw > 1)
    % START SOURCE == ROTLETS AND STOKESLETS
    % START TARGETS == FIBRES
    fsr = 0;
    for k = 2:nw % loop over all walls, except outer wall
        fsr = fsr + o.RSlets(geom.X, walls.center(:,k), lambda(:,k-1), xi(k-1));
    end 
    
    % END TARGETS == FIBRES
    
    % START TARGETS == WALLS
    wsr = zeros(2*Nw,nw);
    for k = 2:nw % loop over all walls, except outer wall
        wsr = wsr + o.RSlets(walls.X, walls.center(:,k), lambda(:,k-1), xi(k-1));
    end   
    
    % END TARGETS == WALLS
    % END SOURCE == ROTLETS
    
    % EVALUATE ROTLET AND STOKESLET EQUATIONS
    valLets = o.letsIntegrals([lambda;xi], etaW, walls);
else
    fsr = 0;
    wsr = 0;
    valLets = [];
end
% EVALUATE VELOCITY ON FIBERS
valFibers = valFibers + ffdlp + fwdlp + fsr;

for k = 1:np
  valFibers(1:end/2,k) = valFibers(1:end/2,k) - Up(1,k);
  valFibers(end/2+1:end,k) = valFibers(end/2+1:end,k) - Up(2,k);
end

for k = 1:np
  valFibers(1:end/2,k) = valFibers(1:end/2,k) ...
                + (geom.X(end/2+1:end,k) - geom.center(2,k))*wp(k);
  valFibers(end/2+1:end,k) = valFibers(end/2+1:end,k)...
                - (geom.X(1:end/2,k) - geom.center(1,k))*wp(k);
end

% EVALUATE VELOCITY ON WALLS
if o.confined
    valWalls = valWalls + wwdlp + wfdlp + wsr;
    valWalls(:,1) = valWalls(:,1)...
                    + pot_walls.exactStokesN0diag(walls, o.N0w, etaW(:,1));
end

% EVALUTATE FORCES ON FIBRES
for k = 1:np
  valForce(1,k) = sum(etaF(1:Np,k).*geom.sa(:,k))*2*pi/Np;
  valForce(2,k) = sum(etaF(Np+1:2*Np,k).*geom.sa(:,k))*2*pi/Np;
end

% EVALUATE TORQUES ON FIBRES
for k = 1:np
  valTorque(k) = sum(((geom.X(Np+1:2*Np,k)-geom.center(2,k)).*etaF(1:Np,k) - ...
                     ((geom.X(1:Np,k)-geom.center(1,k)).*etaF(Np+1:2*Np,k))).*...
                     geom.sa(:,k))*2*pi/Np;
end

% CONSTRUCT OUTPUT VECTOR
Tx = [valFibers(:); valWalls(:); -valForce(:);-valTorque(:);valLets(:)];

if o.profile
    o.om.writeMessage(['Matvec assembly completed in ', ....
                            num2str(toc(tMatvec)), ' seconds']);
end

end % timeMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pz = preconditionerBD(o,z)
% apply the block-diagonal preconditioner whose LU factorization is
% precomputed and stored

Np = (size(o.precop.L,1)-3)/2;
np = size(o.precop.L,3);

if o.confined
    Nw = (size(o.precow.L,1)-3)/2;
    nw = size(o.precow.L,3);
else
    Nw = 0;
    nw = 0;
end

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
% to the stokeslet and rotlet terms.  Center of the rotlet and
% stokeslet is contained in center

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
function z = letsIntegrals(~,otlets,etaM,walls)
% z = letsIntegrals(stokeslet,rotlet,etaM,walls) integrates the density
% function to enforce constraints on stokeslets and rotlets

Nw = walls.N;
nw = walls.n;
z = zeros(3*(nw-1),1);

for k = 2:nw
  stokeslet = otlets(3*(k-2)+1:3*(k-2)+2);
  % two stokeslet terms per inner boundary
  rotlet = otlets(3*(k-1));
  % one rotlet term per inner boundary
  ind = 3*(k-2)+1;
  z(ind) = -2*pi*stokeslet(1) + ...
    sum(etaM(1:Nw,k).*walls.sa(:,k))*2*pi/Nw;
  % integral of density function dotted with [1;0]
  % is one stokeslet
  z(ind+1) = -2*pi*stokeslet(2) + ...
    sum(etaM(Nw+1:2*Nw,k).*walls.sa(:,k))*2*pi/Nw;
  % integral of density fuction dotted with [0;1]
  % is the other stokeslet
  z(ind+2) = -2*pi*rotlet + sum(...
    ((walls.X(Nw+1:2*Nw,k)).*etaM(1:Nw,k) - ...
    (walls.X(1:Nw,k)).*etaM(Nw+1:2*Nw,k)).*...
    walls.sa(:,k))*2*pi/Nw;
  % integral of density function dotted with (-y,x)
  % is the rotlet
end % k

end % letsIntegrals

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
            vInf = [vInf, [-options.couette_speed*y(:,2); options.couette_speed*x(:,2)]];
            
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
            vInf = [-X(end/2+1:end,:);zeros(Np,np)];

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
    
    [Np,np] = size(geom.X);
    [Nw,nw] = size(walls.X);
    
    Ntotal = Np*np + Nw*nw + 3*np + 3*(nw-1);
    
    I = eye(Ntotal);
    M = zeros(Ntotal);
    
    for k = 1:Ntotal
       M(:,k) = o.timeMatVec(I(:,k),geom,walls);
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


