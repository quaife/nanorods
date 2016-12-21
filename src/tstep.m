classdef tstep < handle
% This class defines the functions required to advance the geometry
% forward in time.  Routines that we may need to add are different
% integrators such as semi-implicit, SDC, Runge-Kutta, adaptivity

properties

order        % time stepping order
dt           % time step size
Df           % Stokes double-layer potential for fiber-fiber interaction
Dw           % Stokes double-layer potential for wall-wall interaction
Dupf         % Upsampled Stokes double-layer potential matrix for fibers
Dupw         % Upsampled Stokes double-layer potential matrix for walls
N0w          % N0 matrix to remove rank 1 nullspace
rhs          % Right hand side
ifmm         % flag for using the FMM
inear        % flag for using near-singular integration
gmresTol     % GMRES tolerance
nearStructff % near-singular integration structure (fibre-fibre)
nearStructfw % near-singular integration structure (fibre-wall)
nearStructwf % near-singular integration structure (wall-fibre)
farField     % background flow
usePreco     % use a block-diagonal preconditioner
precoF       % block-diagonal preconditioner for fibres
precoW       % block-diagonal preconditioner for walls
opFibers     % class for fiber layer potentials
opWall       % class for wall layer potentials
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

o.order = options.tstep_order;        
o.ifmm = options.ifmm;                  
o.inear = options.inear;                
o.usePreco = options.usePreco;           
o.profile = options.profile;             
o.confined = options.confined;

o.dt = prams.T/prams.m;                  
o.gmresTol = prams.gmresTol;             

o.farField = @(X) o.bgFlow(X,options); 
o.om = om;
o.tau0 = tau0;

N = prams.N;
Nbd = prams.Nbd;
nv = prams.nv;
nbd = prams.nbd;

% CREATE CLASSES TO EVALUATE POTENTIALS ON FIBRES AND WALLS
o.opFibers = poten(prams.N, om);

if options.confined
    o.opWall = poten(prams.Nbd, om);
    o.Dw = o.opWall.stokesDLmatrix(walls);
else
    o.opWall = [];
    o.Dw = [];
end

% CREATE MATRICES FOR FIBRE-FIBRE SELF INTERATIONS AND 
% WALL-WALL SELF INTERACTIONS
o.Df = o.opFibers.stokesDLmatrix(geom);

% CREATE UPSAMPLED MATRICES
% FIBRE-FIBRE
Xsou = geom.X; 
Nup = N*ceil(sqrt(N));

Xup = [interpft(Xsou(1:N,:),Nup);...
       interpft(Xsou(N+1:2*N,:),Nup)];

geomUp = capsules([],Xup);
o.Dupf = o.opFibers.stokesDLmatrix(geomUp);

% WALL-WALL
if options.confined
    Xsou = walls.X; 
    Nup = Nbd*ceil(sqrt(Nbd));

    Xup = [interpft(Xsou(1:Nbd,:),Nup);...
           interpft(Xsou(Nbd+1:2*Nbd,:),Nup)];

    wallsUp = capsules([],Xup);
    o.Dupw = o.opFibers.stokesDLmatrix(wallsUp);

    % CREATE N0 MATRIX
    o.N0w = o.opWall.stokesN0matrix(walls);
end

% CREATE BLOCK-DIAGONAL PRECONDITIONER
if o.usePreco
    
    if o.profile
        tic;
    end
    
    % FIBEE-FIBRE PRECONDITIONER
    o.precoF.L = zeros(2*N+3,2*N+3,nv);
        o.precoF.U = zeros(2*N+3,2*N+3,nv);
        for k = 1:nv
            [o.precoF.L(:,:,k),o.precoF.U(:,:,k)] =...
                lu([-1/2*eye(2*N)+o.Df(:,:,k) ...
                [-ones(N,1);zeros(N,1)] ...
                [zeros(N,1);-ones(N,1)] ...
                [geom.X(end/2+1:end,k)-geom.center(2,k); -geom.X(1:end/2,k)+geom.center(1,k)];...
                [-geom.sa(:,k)'*2*pi/N zeros(1,N) 0 0 0];
                [zeros(1,N) -geom.sa(:,k)'*2*pi/N 0 0 0];
                [(-geom.X(end/2+1:end,k)'+geom.center(2,k)).*geom.sa(:,k)'*2*pi/N ...
                (geom.X(1:end/2,k)'-geom.center(1,k)).*geom.sa(:,k)'*2*pi/N 0 0 0]]);
        end
        
    if o.confined
         
        % WALL-WALL PRECONDITIONER
        o.precoW.L = zeros(2*Nbd + 3,2*Nbd + 3,nbd);
        o.precoW.U = zeros(2*Nbd + 3,2*Nbd + 3,nbd);
        
        oc = curve;
        sa = walls.sa;
        [x,y] = oc.getXY(walls.X);
        [cx,cy] = oc.getXY(walls.center);
        
        for k = 1:nbd
            
            if k == 1 % first wall does not need Rotlets and Stokeslets
                [o.precoW.L(1:2*Nbd,1:2*Nbd,k),o.precoW.U(1:2*Nbd,1:2*Nbd,k)] =...
                    lu(-1/2*eye(2*Nbd)+o.Dw(:,:,k)+o.N0w(:,:,k));
            
            else  
                r = [x(:,k) - cx(k), y(:,k) - cy(k)];
                rho2 = (x(:,k) - cx(k)).^2 + (y(:,k) - cy(k)).^2;                
                
                col_stokes1 = [-0.5*log(rho2) + r(:,1).*r(:,1)./rho2; r(:,2).*r(:,1)./rho2]/(4*pi);
                col_stokes2 = [r(:,2).*r(:,1)./rho2; -0.5*log(rho2) + r(:,2).*r(:,2)./rho2]/(4*pi);
                col_rot = [r(:,2)./rho2; -r(:,1)./rho2];
                
                int_stokes = 2*pi*sa(:,k)'/Nbd;
                int_rot = 2*pi*[(y(:,k).*sa(:,k))', -(x(:,k).*sa(:,k))']/Nbd;
                
                [o.precoW.L(:,:,k),o.precoW.U(:,:,k)] =...
                    lu([-1/2*eye(2*Nbd)+o.Dw(:,:,k), col_stokes1, col_stokes2, col_rot;...
                    int_stokes, zeros(1,Nbd), -2*pi, 0, 0;...
                    zeros(1,Nbd), int_stokes, 0, -2*pi, 0;...
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
    rhs = o.farField(walls.X);
    o.rhs = [zeros(2*N*nv,1); rhs; zeros(3*nv,1);zeros(3*(nbd-1),1)];
else
    rhs = o.farField(geom.X);
    o.rhs = [-rhs(:); zeros(2*Nbd*nbd,1); zeros(3*nv,1)];
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
    o.farField = @(X) o.bgFlow(X,options); 

    N = prams.N;
    Nbd = prams.Nbd;
    nv = prams.nv;
    nbd = prams.nbd;
    
    rhs = o.farField(geom.X);
    o.rhs = [-rhs(:); zeros(2*Nbd*nbd,1); zeros(3*nv,1)];
end

% CREATE NEAR SINGULAR INTEGRATION STRUCTURES
if o.inear
    if o.profile
        tic;
    end
    
    o.nearStructff = geom.getZone([],1);
    
    if o.confined
        [~,o.nearStructfw] = geom.getZone(walls,2);
        [~,o.nearStructwf] = walls.getZone(geom,2);
    end
    
    if o.profile
        o.om.writeMessage(['getZone ... ', num2str(toc)]);
    end
end

% ROTATE FIBRE-FIBRE DLP AND FIBRE-FIBRE PRECONDITIONER
N = geom.N;
nv = geom.nv;

if options.confined
    Nbd = walls.N;
    nbd = walls.nv;
else
    Nbd = 0;
    nbd = 0;
end

Nup = N*ceil(sqrt(N));
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

for i = 1:nv

    R = spdiags([sin(dtau(i))*ones(2*N,1), cos(dtau(i))*ones(2*N,1)...
                -sin(dtau(i))*ones(2*N,1)], [-N, 0, N], zeros(2*N, 2*N));
    
    Rup = spdiags([sin(dtau(i))*ones(2*Nup,1), cos(dtau(i))*ones(2*Nup,1)...
                -sin(dtau(i))*ones(2*Nup,1)], [-Nup, 0, Nup], zeros(2*Nup, 2*Nup));  
            
    o.Df(:,:,i) = R*o.Df(:,:,i)*R';
    o.Dupf(:,:,i) = Rup*o.Dupf(:,:,i)*Rup';
    
    if o.usePreco
        
        [o.precoF.L(:,:,i),o.precoF.U(:,:,i)] =...
            lu([-1/2*eye(2*N)+o.Df(:,:,i) ...
            [-ones(N,1);zeros(N,1)] ...
            [zeros(N,1);-ones(N,1)] ...
            [geom.X(end/2+1:end,i)-geom.center(2,i); -geom.X(1:end/2,i)+geom.center(1,i)];...
            [-geom.sa(:,i)'*2*pi/N zeros(1,N) 0 0 0];
            [zeros(1,N) -geom.sa(:,i)'*2*pi/N 0 0 0];
            [(-geom.X(end/2+1:end,i)'+geom.center(2,i)).*geom.sa(:,i)'*2*pi/N ...
            (geom.X(1:end/2,i)'-geom.center(1,i)).*geom.sa(:,i)'*2*pi/N 0 0 0]]);
    end    
end

% SOLVE SYSTEM USING GMRES
maxit = 2*N*nv;% should be a lot lower than this
%o.rhs  = [ones(2*nv*N,1);2*ones(2*Nbd*nbd,1);3*ones(2*nv,1);4*ones(nv,1);5*ones(2*(nbd-1),1);6*ones(nbd-1,1)];
% 
% clf;
% semilogy(abs(o.timeMatVec(o.preconditionerBD(o.rhs), geom, walls)-o.rhs), '.');

if o.usePreco
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom,walls),o.rhs,[],o.gmresTol,...
      maxit,@o.preconditionerBD,[], o.rhs);
else
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom,walls),o.rhs,[],o.gmresTol,...
      maxit);
end

iter = I(2);

% REORGANIZE COLUMN VECTOR INTO MATRIX
% EXTRACT DENSITY FUNCITONS ON FIBRES AND WALLS
% each column of etaF corresponds to the density function of a rigid body
etaF = zeros(2*N,nv);
for k = 1:nv
  etaF(:,k) = Xn((k-1)*2*N+1:k*2*N);
end

% each column of etaW corresponds to the density function of a solid wall
if options.confined
    etaW = zeros(2*Nbd, nbd);
    for k = 1:nbd
       etaW(:,k) = Xn(2*N+1+(k-1)*2*Nbd:2*N+k*2*Nbd);
    end
else
    etaW = 0;
end

% EXTRACT TRANSLATIONAL AND ROTATIONAL VELOCITIES
Up = zeros(2,nv);
for k = 1:nv
  Up(:,k) = Xn(2*N*nv+2*Nbd*nbd+(k-1)*2+1:2*N*nv+2*Nbd*nbd+k*2);
end

wp = zeros(1,nv);
for k = 1:nv
  wp(k) = Xn(2*N*nv+2*Nbd*nbd+2*nv+k);
end

if options.confined
    stokes = zeros(2,nbd-1);
    for k = 1:nbd-1
       stokes(:,k) = Xn(2*N*nv+2*Nbd*nbd+3*nv+1+2*(k-1): 2*N*nv+2*Nbd*nbd+3*nv+2*k);
    end

    rot = zeros(1,nbd-1);
    for k = 1:nbd-1
       rot(k) = Xn(2*N*nv+2*Nbd*nbd+3*nv+2*(nbd-1)+k);
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

N = geom.N;   % points per body
nv = geom.nv; % number of bodies

if o.confined
    Nbd = walls.N;   % points per wall
    nbd = walls.nv; % number of wallls
else
    Nbd = 0;
    nbd = 0;
end

potFibers = o.opFibers;
potWalls = o.opWall;

% Output of Tx that corresponds to the velocity of the fibers
valFibers = zeros(2*N,nv);
% Output of Tx that corresponds to the velocity of the walls
valWalls = zeros(2*Nbd,nbd);
% output of Tx that corresponds to force on fibres
valForce = zeros(2,nv);
% output of Tx that corresponds to torque on fibres
valTorque = zeros(nv,1);

% BEGIN FORMATTING UNKNOWN VECTOR

% EXTRACT DENSITY FUNCTION FOR FIBRES
etaF = zeros(2*N,nv);
for k = 1:nv
  etaF(:,k) = Xn((k-1)*2*N+1:k*2*N);
end

% EXTRACT DENSITY FUNCTION FOR WALLS
etaW = zeros(2*Nbd,nbd);
if o.confined
    for k = 1:nbd
      etaW(:,k) = Xn(2*N*nv+1+(k-1)*2*Nbd:2*N*nv+k*2*Nbd);
    end
end

% EXTRACT TRANSLATIONAL AND ROTATIONAL VELOCITIES OF FIBRES
Up = zeros(2,nv);
for k = 1:nv
  Up(:,k) = Xn(2*N*nv+2*Nbd*nbd+1+2*(k-1):2*N*nv+2*Nbd*nbd+2*k);
end
wp = zeros(1,nv);
for k = 1:nv
  wp(k) = Xn(2*N*nv+2*Nbd*nbd+2*nv+k);
end

if o.confined
    % EXTRACT ROTLETS AND STOKESLETS
    lambda = zeros(2,nbd-1);
    for k = 1:nbd-1
       lambda(:,k) = Xn(2*N*nv+2*Nbd*nbd+3*nv+1+2*(k-1):2*N*nv+2*Nbd*nbd+3*nv+2*k); 
    end

    xi=zeros(1,nbd-1);
    for k=1:nbd-1
        xi(k) = Xn(2*N*nv+2*Nbd*nbd+3*nv+2*(nbd-1)+k);
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
valFibers = valFibers + potFibers.exactStokesDLdiag(geom, o.Df, etaF);
if o.confined
    valWalls = valWalls + potWalls.exactStokesDLdiag(walls, o.Dw, etaW);
end

% START OF SOURCE == FIBRES
% START OF TARGET == FIBRES
if o.ifmm
  kernel = @potFibers.exactStokesDLfmm;
else
  kernel = @potFibers.exactStokesDL;
end

kernelDirect = @potFibers.exactStokesDL;

if o.inear
  DLP = @(X) potFibers.exactStokesDLdiag(geom,o.Df,X) - 1/2*X;
  ffdlp = potFibers.nearSingInt(geom, etaF, DLP, o.Dupf, o.nearStructff ,kernel,...
        kernelDirect, geom, true, false);
else
  ffdlp = kernel(geom, etaF, o.Df);
end
% END OF TARGET == FIBRES

% START OF TARGET == WALLS 
if o.confined
   
   if o.ifmm
       kernel = @potWalls.exactStokesDLfmm;       
   else
       kernel = @potWalls.exactStokesDL;
   end
   
   kernelDirect = @potWalls.exactStokesDL;
   
   if o.inear
       DLP = @(X) potWalls.exactStokesDLdiag(geom, o.Df, X) - 1/2*X;
       wfdlp = potWalls.nearSingInt(geom, etaF, DLP, o.Dupf, o.nearStructfw, ...
           kernel, kernelDirect, walls, false, false);
   else
       wfdlp = kernel(geom, etaF);
   end
else
    wfdlp = zeros(2*Nbd,nbd);
end
% END OF TARGET == WALLS

% START OF SOURCE == WALLS
% START OF TARGET == FIBERS
if o.confined
   
   if o.ifmm
        kernel = @potFibers.exactStokesDLfmm;       
   else
       kernel = @potFibers.exactStokesDL;
   end
   
   kernelDirect = @potFibers.exactStokesDL;
   
   if o.inear
       DLP = @(X) potFibers.exactStokesDLdiag(walls, o.Dw, X) - 1/2*X;
       fwdlp = potFibers.nearSingInt(walls, etaW, DLP, o.Dupw, o.nearStructwf, ...
           kernel, kernelDirect, geom, false, false);
   else
       fwdlp = kernel(walls, etaW);
   end
else
    fwdlp = zeros(2*N,nv);
end
% END OF TARGET == FIBRES

% START OF TARGET == WALLS
if o.confined
    
    if o.ifmm
        kernel = @potWalls.exactStokesDLfmm;
    else
        kernel = @potWalls.exactStokesDL;
    end
    
    wwdlp = kernel(walls, etaW, o.Dw);
    
else
    wwdlp = zeros(2*Nbd,nbd);
end
% END OF TARGET == WALLS
% END OF SOURCE == WALLS

if (o.confined && nbd > 1)
    % START SOURCE == ROTLETS AND STOKESLETS
    % START TARGETS == FIBRES
    fsr = 0;
    for k = 2:nbd % loop over all walls, except outer wall
        fsr = fsr + o.RSlets(geom.X, walls.center(:,k), lambda(:,k-1), xi(k-1));
    end 
    
    % END TARGETS == FIBRES
    
    % START TARGETS == WALLS
    wsr = zeros(2*Nbd,nbd);
    for k = 2:nbd % loop over all walls, except outer wall
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

for k = 1:nv
  valFibers(1:end/2,k) = valFibers(1:end/2,k) - Up(1,k);
  valFibers(end/2+1:end,k) = valFibers(end/2+1:end,k) - Up(2,k);
end

for k = 1:nv
  valFibers(1:end/2,k) = valFibers(1:end/2,k) + (geom.X(end/2+1:end,k) - geom.center(2,k))*wp(k);
  valFibers(end/2+1:end,k) = valFibers(end/2+1:end,k) - (geom.X(1:end/2,k) - geom.center(1,k))*wp(k);
end

% EVALUATE VELOCITY ON WALLS
if o.confined
    valWalls = valWalls + wwdlp + wfdlp + wsr;
    valWalls(:,1) = valWalls(:,1) + potWalls.exactStokesN0diag(walls, o.N0w, etaW(:,1));
end

% EVALUTATE FORCES ON FIBRES
for k = 1:nv
  valForce(1,k) = sum(etaF(1:N,k).*geom.sa(:,k))*2*pi/N;
  valForce(2,k) = sum(etaF(N+1:2*N,k).*geom.sa(:,k))*2*pi/N;
end

% EVALUATE TORQUES ON FIBRES
for k = 1:nv
  valTorque(k) = sum(((geom.X(N+1:2*N,k)-geom.center(2,k)).*etaF(1:N,k) - ...
                     ((geom.X(1:N,k)-geom.center(1,k)).*etaF(N+1:2*N,k))).*...
                     geom.sa(:,k))*2*pi/N;
end

% CONSTRUCT OUTPUT VECTOR
Tx = [valFibers(:); valWalls(:); -valForce(:);-valTorque(:);valLets(:)];

if o.profile
    o.om.writeMessage(['Matvec assembly completed in ', num2str(toc(tMatvec)), ' seconds']);
end

end % timeMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pz = preconditionerBD(o,z)
% apply the block-diagonal preconditioner whose LU factorization is
% precomputed and stored

N = (size(o.precoF.L,1)-3)/2;
nv = size(o.precoF.L,3);

if o.confined
    Nbd = (size(o.precoW.L,1)-3)/2;
    nbd = size(o.precoW.L,3);
else
    Nbd = 0;
    nbd = 0;
end

Pz = z;

% APPLY FIBRE-FIBRE PRECONDITIONER
for k = 1:nv
  etaStart = (k-1)*2*N + 1;
  endEnd = etaStart + 2*N - 1;
  uStart = 2*N*nv + 2*Nbd*nbd + (k-1)*2 + 1;
  uEnd = uStart + 1;
  omegaStart = 2*N*nv + 2*Nbd*nbd + nv*2 + (k-1) + 1;
  omegaEnd = omegaStart;
  zblock = o.precoF.U(:,:,k)\(o.precoF.L(:,:,k)\...
      [z(etaStart:endEnd);z(uStart:uEnd);z(omegaStart:omegaEnd)]);

  Pz(etaStart:endEnd) = zblock(1:2*N);
  Pz(uStart:uEnd) = zblock(2*N+1:2*N+2);
  Pz(omegaStart:omegaEnd) = zblock(2*N+3:2*N+3);
end

% APPLY WALL-WALL PRECONDITIONER
for k = 1:nbd
    
    if k == 1 %No stokeslets or rotlets
       xiStart = 2*N*nv +  1;
       xiEnd = xiStart + 2*Nbd - 1;
       
       zblock = o.precoW.U(1:2*Nbd,1:2*Nbd,1)\(o.precoW.L(1:2*Nbd,1:2*Nbd,1)\...
                            z(xiStart:xiEnd));
       Pz(xiStart:xiEnd) = zblock;
       
    else
        xiStart = 2*N*nv + 2*(k-1)*Nbd + 1;
        xiEnd = xiStart + 2*Nbd - 1;
        
        stokesletStart = 2*N*nv+2*Nbd*nbd+3*nv+2*(k-2)+1;
        stokesletEnd = stokesletStart + 1;
        rotletStart = 2*N*nv+2*Nbd*nbd+3*nv+2*(nbd-1)+(k-1);
        rotletEnd = rotletStart;
        
        zblock = o.precoW.U(:,:,k)\(o.precoW.L(:,:,k)\...
                [z(xiStart:xiEnd);z(stokesletStart:stokesletEnd);...
                z(rotletStart:rotletEnd)]);

        Pz(xiStart:xiEnd) = zblock(1:2*Nbd);
        Pz(stokesletStart:stokesletEnd) = zblock(2*Nbd+1:2*Nbd+2);
        Pz(rotletStart:rotletEnd) = zblock(2*Nbd+3:2*Nbd+3);
    end
end

end % preconditioner

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = RSlets(o,X,center,stokeslet,rotlet)
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
function z = letsIntegrals(o,otlets,etaM,walls)
% z = letsIntegrals(stokeslet,rotlet,etaM,walls) integrates the density
% function to enforce constraints on stokeslets and rotlets

Nbd = walls.N;
nbd = walls.nv;
z = zeros(3*(nbd-1),1);

for k = 2:nbd
  stokeslet = otlets(3*(k-2)+1:3*(k-2)+2);
  % two stokeslet terms per inner boundary
  rotlet = otlets(3*(k-1));
  % one rotlet term per inner boundary
  ind = 3*(k-2)+1;
  z(ind) = -2*pi*stokeslet(1) + ...
    sum(etaM(1:Nbd,k).*walls.sa(:,k))*2*pi/Nbd;
  % integral of density function dotted with [1;0]
  % is one stokeslet
  z(ind+1) = -2*pi*stokeslet(2) + ...
    sum(etaM(Nbd+1:2*Nbd,k).*walls.sa(:,k))*2*pi/Nbd;
  % integral of density fuction dotted with [0;1]
  % is the other stokeslet
  z(ind+2) = -2*pi*rotlet + sum(...
    ((walls.X(Nbd+1:2*Nbd,k)).*etaM(1:Nbd,k) - ...
    (walls.X(1:Nbd,k)).*etaM(Nbd+1:2*Nbd,k)).*...
    walls.sa(:,k))*2*pi/Nbd;
  % integral of density function dotted with (-y,x)
  % is the rotlet
end % k

end % letsIntegrals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vInf = bgFlow(o,X,options)
    
N = size(X,1)/2;
nv = size(X,2);
oc = curve;

[x,y] = oc.getXY(X);

if options.confined
    switch options.farField
        case 'constant'
            vInf = [ones(N,nv);zeros(N,nv)];

        case 'circle'
            vInf = [-y(:,1); x(:,1)];
            
        case 'couette'
            %[-y,x]==>clockwise
            %vInf = 1*[-y(:,1)+mean(y(:,1));x(:,1)-mean(x(:,1))];
            vInf = zeros(2*N,1);
            %[y,-x]==>counter clockwise
            vInf = [vInf; [-y(:,2)+mean(y(:,2)); x(:,2)-mean(x(:,2))]];

        otherwise
            vInf = [zeros(N,nv);zeros(N,nv)];
    end
else
    switch options.farField
        case 'shear'
            vInf = [5*X(end/2+1:end,:);zeros(N,nv)];

        case 'extensional'
            vInf = [X(1:end/2,:);-X(end/2+1:end,:)];

        otherwise
            vInf = [ones(N,nv);zeros(N,nv)];
    end

end
    
end % bgFlow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = build_matrix(o, geom, walls)
    
    [N,nv] = size(geom.X);
    [Nbd,nbd] = size(walls.X);
    
    Ntotal = N*nv + Nbd*nbd + 3*nv + 3*(nbd-1);
    
    I = eye(Ntotal);
    M = zeros(Ntotal);
    
    for k = 1:Ntotal
       M(:,k) = o.timeMatVec(I(:,k),geom,walls);
    end    
end % build_matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = build_preconditioner_matrix(o, geom, walls)

    [N,nv] = size(geom.X);
    [Nbd,nbd] = size(walls.X);
    
    N = N/2;
    Nbd = Nbd/2;
    
    Ntotal = 2*N*nv + 2*Nbd*nbd + 3*nv + 3*(nbd-1);
    
    M = zeros(Ntotal,Ntotal);
    
    % fibre-fibre
    for k = 1:nv
        start_row = 1+(k-1)*2*N;
       M(start_row:start_row+2*N-1,start_row:start_row+2*N-1) =  -1/2*eye(2*N)+o.Df(:,:,k);
    end
    
    % wall-wall
    for k = 1:nbd
       start_row = 2*N*nv+1+(k-1)*2*Nbd;
       M(start_row:start_row+2*Nbd-1,start_row:start_row+2*Nbd-1) = -1/2*eye(2*Nbd)+o.Dw(:,:,k)+o.N0w(:,:,k);    
    end
    
    % u and omega    
    for k = 1:nv
        
       start_row = 2*N*nv+2*Nbd*nbd+2*(k-1)+1;
       start_col = 2*N*(k-1)+1;
       
       M(start_row:start_row+1,start_col:start_col+2*N-1) = ...
                [[-geom.sa(:,k)'*2*pi/N, zeros(1,N)];...
                [zeros(1,N), -geom.sa(:,k)'*2*pi/N]];
       
       start_row = 2*N*nv+2*Nbd*nbd+2*nv+k;
       M(start_row,start_col:start_col+2*N-1) = ...
                 [(-geom.X(end/2+1:end,k)'+geom.center(2,k)).*geom.sa(:,k)'*2*pi/N ...
                    (geom.X(1:end/2,k)'-geom.center(1,k)).*geom.sa(:,k)'*2*pi/N];
        
       start_row = 2*N*(k-1)+1;
       start_col = 2*N*nv+2*Nbd*nbd+2*(k-1)+1;
       
       M(start_row:start_row+2*N-1,start_col:start_col+1) = ...
                    [[-ones(N,1);zeros(N,1)] ...
                    [zeros(N,1);-ones(N,1)]];
       
        start_col = 2*N*nv + 2*Nbd*nbd+2*nv+k;
        M(start_row:start_row+2*N-1,start_col) = ...
                [geom.X(end/2+1:end,k)-geom.center(2,k); -geom.X(1:end/2,k)+geom.center(1,k)];
    end
    
    
    % stokeslets and rotlets
    oc = curve;
    sa = walls.sa;
    [x,y] = oc.getXY(walls.X);
    [cx,cy] = oc.getXY(walls.center);
    
    for k = 1:nbd-1
        r = [x(:,k+1) - cx(k+1), y(:,k+1) - cy(k+1)];
        rho2 = (x(:,k+1) - cx(k+1)).^2 + (y(:,k+1) - cy(k+1)).^2;
        
        col_stokes1 = [-0.5*log(rho2) + r(:,1).*r(:,1)./rho2; r(:,2).*r(:,1)./rho2]/(4*pi);
        col_stokes2 = [r(:,2).*r(:,1)./rho2; -0.5*log(rho2) + r(:,2).*r(:,2)./rho2]/(4*pi);
        col_rot = [r(:,2)./rho2; -r(:,1)./rho2];
        
        int_stokes = 2*pi*sa(:,k+1)'/Nbd;
        int_rot = 2*pi*[(y(:,k+1).*sa(:,k+1))', -(x(:,k+1).*sa(:,k+1))']/Nbd;
        
        start_row = 2*N*nv+2*Nbd*nbd+3*nv+1;
        start_col = 2*N*nv+2*Nbd*k+1;
        
        M(start_row:start_row+2,start_col:start_col+2*Nbd-1) = ...
                    [int_stokes, zeros(1,Nbd);...
                    zeros(1,Nbd), int_stokes;...
                    int_rot];
        
        start_col = 2*N*nv+2*Nbd*nbd+3*nv+(k-1)+1;
        M(start_row:start_row+2,start_col:start_col+2) = -2*pi*eye(3);
        
        start_row =  2*N*nv+2*Nbd*k+1;
        start_col = 2*N*nv+2*Nbd*nbd+3*nv+1 +3*(k-1);
        
        M(start_row:start_row+2*Nbd-1,start_col:start_col+2) = [col_stokes1, col_stokes2, col_rot];
    end
    
end % build_preconditioner_matrix

end % methods

end % classdef


