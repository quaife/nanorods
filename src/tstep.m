classdef tstep < handle
% This class defines the functions required to advance the geometry
% forward in time.  Routines that we may need to add are different
% integrators such as semi-implicit, SDC, Runge-Kutta, adaptivity

properties

order      % time stepping order
dt         % time step size
Df         % Stokes double-layer potential for fiber-fiber interaction
Dw         % Stokes double-layer potential for wall-wall interaction
Dupf       % Upsampled Stokes double-layer potential matrix for fibers
Dupw       % Upsampled Stokes double-layer potential matrix for walls
rhs        % Right hand side
ifmm       % flag for using the FMM
inear      % flag for using near-singular integration
gmresTol   % GMRES tolerance
nearStruct % near-singular integration structure 
farField   % background flow
usePreco   % use a block-diagonal preconditioner
bdiagPreco % block-diagonal preconditioner
opFibers   % class for fiber layer potentials
opWall     % class for wall layer potentials
profile    % flag to time certain parts of code
om         % monitor class

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(options, prams, om, geom, walls)
% o.tstep(options,prams): constructor.  Initialize class.  Take all
% elements of options and prams needed by the time stepper

o.order = options.tstep_order;        
o.ifmm = options.ifmm;                  
o.inear = options.inear;                
o.usePreco = options.usePreco;           
o.profile = options.profile;             

o.dt = prams.T/prams.m;                  
o.gmresTol = prams.gmresTol;             

o.farField = @(X) o.bgFlow(X,options); 
o.om = om;

N = prams.N;
Nbd = prams.Nbd;

% CREATE CLASSES TO EVALUATE POTENTIALS ON FIBRES AND WALLS
o.opFibers = poten(prams.N, om);

if options.confined
    o.opWall = poten(prams.Nbd, om);
else
    o.opWall = [];
end

% CREATE MATRICES FOR FIBRE-FIBRE SELF INTERATIONS AND 
% WALL-WALL SELF INTERACTIONS
o.Df = o.opFibers.stokesDLmatrix(geom);
o.Dw = o.opWall.stokesDLmatrix(walls);

% CREATE UPSAMPLES MATRICES
% FIBRE-FIBRE
Xsou = geom.X; 
Nup = N*ceil(sqrt(N));

Xup = [interpft(Xsou(1:N,:),Nup);...
       interpft(Xsou(N+1:2*N,:),Nup)];

geomUp = capsules([],Xup);
o.Dupf = o.opFibers.stokesDLmatrix(geomUp);

% WALL-WALL
Xsou = walls.X; 
Nup = Nsou*ceil(sqrt(Nbd));

Xup = [interpft(Xsou(1:Nbd,:),Nup);...
       interpft(Xsou(Nbd+1:2*Nbd,:),Nup)];

wallsUp = capsules([],Xup);
o.Dupw = o.opFibers.stokesDLmatrix(wallsUp);

% CREATE BLOCK-DIAGONAL PRECONDITIONER
if o.usePreco
    
    if o.profile
        tic;
    end
    
    if o.confined
        o.bdiagPreco.L = zeros(2*N+3,2*N+3,nv);
        o.bdiagPreco.U = zeros(2*N+3,2*N+3,nv);
        for k = 1:nv
            [o.bdiagPreco.L(:,:,k),o.bdiagPreco.U(:,:,k)] =...
                lu([-1/2*eye(2*N)+o.D(:,:,k) ...
                [-ones(N,1);zeros(N,1)] ...
                [zeros(N,1);-ones(N,1)] ...
                [geom.X(end/2+1:end,k)-geom.center(2,k); -geom.X(1:end/2,k)+geom.center(1,k)];...
                [-geom.sa(:,k)'*2*pi/N zeros(1,N) 0 0 0];
                [zeros(1,N) -geom.sa(:,k)'*2*pi/N 0 0 0];
                [(-geom.X(end/2+1:end,k)'+geom.center(2,k)).*geom.sa(:,k)'*2*pi/N ...
                (geom.X(1:end/2,k)'-geom.center(1,k)).*geom.sa(:,k)'*2*pi/N 0 0 0]]);
            
            % whole diagonal block which doesn't have a null space
        end
        
    else
        o.bdiagPreco.L = zeros(2*N+3,2*N+3,nv);
        o.bdiagPreco.U = zeros(2*N+3,2*N+3,nv);
        for k = 1:nv
            [o.bdiagPreco.L(:,:,k),o.bdiagPreco.U(:,:,k)] =...
                lu([-1/2*eye(2*N)+o.D(:,:,k) ...
                [-ones(N,1);zeros(N,1)] ...
                [zeros(N,1);-ones(N,1)] ...
                [geom.X(end/2+1:end,k)-geom.center(2,k); -geom.X(1:end/2,k)+geom.center(1,k)];...
                [-geom.sa(:,k)'*2*pi/N zeros(1,N) 0 0 0];
                [zeros(1,N) -geom.sa(:,k)'*2*pi/N 0 0 0];
                [(-geom.X(end/2+1:end,k)'+geom.center(2,k)).*geom.sa(:,k)'*2*pi/N ...
                (geom.X(1:end/2,k)'-geom.center(1,k)).*geom.sa(:,k)'*2*pi/N 0 0 0]]);
            
            % whole diagonal block which doesn't have a null space
        end
    end
    
    if o.profile
        o.om.writeMessage(['Building preconditioner ... ', num2str(toc)]);
    end
    
    % CREATE RIGHT HAND SIDE
if options.confined
    rhs = o.farField(walls.X);
    o.rhs = [zeros(2*N*nv,1); rhs; zeros(3*nv,1)];    
else
   rhs = o.farField(geom.X); 
   o.rhs = [-rhs; zeros(2*Nbd*nvbd,1); zeros(3*nv,1)];
end
  
end % constructor: tstep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eta,Up,wp,iter,iflag,res] = timeStep(o,geom,tau)
% [X,iter,iflag] = timeStep(Xstore) takes the current configuration
% in Xstore (can be a three-dimensional array corresponding to  previous
% time steps if doing multistep) and returns a new shape X, the number
% of GMRES iterations, and the GMRES flag 

if o.inear
    if o.profile
        tic;
    end
    
    o.nearStruct = geom.getZone([],1);
    
    if o.profile
        o.om.writeMessage(['getZone ... ', num2str(toc)]);
    end
end

% ROTATE FIBRE-FIBRE DLP AND FIBRE-FIBRE PRECONDITIONER
N = geom.N;
nv = geom.nv;

Nup = N*ceil(sqrt(N));

for i = 1:nv
    R = spdiags([sin(tau(i))*ones(2*N,1), cos(tau(i))*ones(2*N,1)...
                -sin(tau(i))*ones(2*N,1)], [-N, 0, N], zeros(2*N, 2*N));
    
    Rup = spdiags([sin(tau(i))*ones(2*Nup,1), cos(tau(i))*ones(2*Nup,1)...
                -sin(tau(i))*ones(2*Nup,1)], [-Nup, 0, Nup], zeros(2*Nup, 2*Nup));  
            
    o.Df(:,:,i) = R*o.Df(:,:,i)*R';
    o.bdiagPreco.L(1:2*N,1:2*N,i) = R*o.bdiagPreco.L(1:2*N,1:2*N,i)*R';
    o.bdiagPreco.U(1:2*N,1:2*N,i) = R*o.bdiagPreco.U(1:2*N,1:2*N,i)*R';
    
    o.Dupf(:,:,i) = Rup*o.Dupf(:,:,i)*Rup';
end

% SOLVE SYSTEM USING GMRES
maxit = 2*N*nv;% should be a lot lower than this
if ~o.usePreco
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom),o.rhs,[],o.gmresTol,...
      maxit);
else
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom),o.rhs,[],o.gmresTol,...
      maxit,@o.preconditionerBD);
end

iter = I(2);

% REORGANIZE COLUMN VECTOR INTO MATRIX
% each column corresponds to the density function of the rigid body
eta = zeros(2*N,nv);
for k = 1:nv
  eta(1:2*N,k) = Xn((k-1)*2*N+1:k*2*N);
end

% EXTRACT TRANSLATIONAL AND ROTATIONAL VELOCITIES
Up = zeros(2,nv);
for k = 1:nv
  Up(:,k) = Xn(2*N*nv+(k-1)*2+1:2*N*nv+k*2);
end

wp = zeros(1,nv);
for k = 1:nv
  wp(k) = Xn(2*N*nv+2*nv+k);
end

end % timeStep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tx = timeMatVec(o,Xn,geom,walls)
% Tx = timeMatVec(Xn,geom) does a matvec for GMRES 
% Xn is a state vector that contains the following in order:
% fiber1: eta_x, eta_y fiber2: eta_x, eta_y ... fiberNv: eta_x, eta_y
% wall1 : xi_x, xi_y wall2: xi_x, xi_y ... wallNvbd: xi_x, xi_y
% fiber1: u, v fiber2: u, v ... fiberNv: u, v
% fiber1: omega fiber2: omega ... fiberNv: omega

if o.profile
    tMatvec = tic;
end

N = geom.N;   % points per body
nv = geom.nv; % number of bodies

if o.confined
    Nbd = walls.N;   % points per wall
    nvbd = walls.nv; % number of wallls
else
    Nbd = 0;
    nv = 0;
end

potFibers = o.opFibers;
potWalls = o.opWall;

% Output of Tx that corresponds to the velocity of the fibers
valFibers = zeros(2*N,nv);
% Output of Tx that corresponds to the velocity of the walls
valWalls = zeros(2*Nbd,nvbd);
% output of Tx that corresponds to force on fibres
valForce = zeros(2,nv);
% output of Tx that corresponds to torque on fibres
valTorque = zeros(nv,1);

% EXTRACT DENSITY FUNCTION FOR FIBRES
etaFibers = zeros(2*N,nv);
for k = 1:nv
  etaFibers(:,k) = Xn((k-1)*2*N+1:k*2*N);
end

% EXTRACT DENSITY FUNCTION FOR WALLS
etaWalls = zeros(2*Nbd,nvbd);
for k = 1:nvbd
  etaWalls(:,k) = Xn(2*N*nv+1+(k-1)*Nbd:2*N*nv+1+k*Nbd);
end

% EXTRACT TRANSLATIONAL AND ROTATIONAL VELOCITIES OF FIBRES
Up = zeros(2,nv);
for k = 1:nv
  Up(:,k) = Xn(2*N*nv+2*Nbd*nvbd+1+2*(k-1):2*N*nv+2*Nbd*nvbd+1+2*k);
end
wp = zeros(1,nv);
for k = 1:nv
  wp(k) = Xn(2*N*nv+2*Nbd*nvbd+2*nv+k);
end

% ADD JUMP IN DLP
valFibers = valFibers - 1/2*etaFibers;
valWalls = valWalls - 1/2*etaWalls;

% ADD SELF CONTRIBUTION
valFibers = valFibers + potFibers.exactStokesDLdiag(geom, o.Df, etaFibers);
valWalls = valWalls + potWalls.exactStokesDLdiag(walls, o.Dw, etaWalls);

% START OF SOURCE == FIBRES

% START OF TARGET == FIBRES
% COMPUTE FIBRE-FIBRE DLP
if o.ifmm
  kernel = @potFibers.exactStokesDLfmm;
else
  kernel = @potFibers.exactStokesDL;
end

kernelDirect = @potFibers.exactStokesDL;

if o.inear
  DLP = @(X) potFibers.exactStokesDLdiag(geom,o.Df,X) - 1/2*X;
  ffdlp = potFibers.nearSingInt(geom, etaFibres, DLP, o.nearStruct ,kernel,...
        kernelDirect, geom, true, false);
else
  ffdlp = kernel(geom, etaFibers);
end
% END OF TARGET == FIBRES

% START OF TARGET == WALLS 
% COMPUTE WALL-FIBRE DLP
if o.confined
   
   if o.ifmm
        kernel = @potFibres.exactStokesDLfmm;       
   else
       kernel = @potFibres.exactStokesDL;
   end
   
   kernelDirect = @potFibres.exactStokesDL;
   
   if o.inear
       DLP = @(X) potFibres.exactStokesDLdiag(geom, o.Df, X);
       wfdlp = potWalls.nearSingInt(walls, etaWalls, DLP, o.nearStruct, ...
           kernel, kernelDirect, geom, false, false);
   else
       wfdlp = kernal(geom, etaWalls);
   end
else
    wfdlp = zeros(2*N,nv);
end
% END OF TARGET == WALLS

% START OF SOURCE == WALLS

% START OF TARGET == FIBERS
% COMPUTE FIBRE-WALL DLP
if o.confined
   
   if o.ifmm
        kernel = @potWall.exactStokesDLfmm;       
   else
       kernel = @potwall.exactStokesDL;
   end
   
   kernelDirect = @potWall.exactStokesDL;
   
   if o.inear
       DLP = @(X) potWalls.exactStokesDLdiag(walls, o.Dw, X);
       fwdlp = potWalls.nearSingInt(geom, etaFibres, DLP, o.nearStruct, ...
           kernel, kernelDirect, walls, false, false);
   else
       fwdlp = kernal(walls, etaFibers);
   end
else
    fwdlp = zeros(2*N,nv);
end

% WALL-WALL INTERACTIONS ARE COMPUTED IN CONSTRUCTOR

% END OF SOURCE == WALLS


% EVALUATE VELOCITY ON FIBERS
valFibers = valFibers + ffdlp + potWalls.exactStokesDL(geom,etaFibers,D,Xtar,K1);

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
    valWalls = valWalls - 1/2*etaWalls + potWalls.exactStokesDLdiag(walls, o.DWall, etaWalls);
    valWalls(:,1) = valWalls(:,1) + potWalls.exactStokesN0diag(walls, o.N0wall, etaWalls(:,1));
    
    valWalls = valWalls + fwdlp;
end

% EVALUTATE FORCES ON FIBRES
for k = 1:nv
  valForce(1,k) = sum(eta(1:N,k).*geom.sa(:,k))*2*pi/N;
  valForce(2,k) = sum(eta(N+1:2*N,k).*geom.sa(:,k))*2*pi/N;
end

% EVALUATE TORQUES ON FIBRES
for k = 1:nv
  valTorque(k) = sum(((geom.X(N+1:2*N,k)-geom.center(2,k)).*eta(1:N,k) - ...
                     ((geom.X(1:N,k)-geom.center(1,k)).*eta(N+1:2*N,k))).*...
                     geom.sa(:,k))*2*pi/N;
end

% CONSTRUCT OUTPUT VECTOR
Tx = [valFibers(:); valWalls(:); -valForce(:);-valTorque(:)];

if o.profile
    o.om.writeMessage(['Matvec assembly completed in ', num2str(toc(tMatvec)), ' seconds']);
end

end % timeMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pz = preconditionerBD(o,z)
% apply the block-diagonal preconditioner whose LU factorization is
% precomputed and stored

N = (size(o.bdiagPreco.L,1)-3)/2;
nv = size(o.bdiagPreco.L,3);

Pz = zeros(2*N*nv+3*nv,1);
for k = 1:nv
  istart1 = (k-1)*2*N + 1;
  iend1 = istart1 + 2*N - 1;
  istart2 = 2*N*nv + (k-1)*2 + 1;
  iend2 = istart2 + 1;
  istart3 = 2*N*nv + nv*2 + (k-1) + 1;
  iend3 = istart3;
  zblock = o.bdiagPreco.U(:,:,k)\...
      (o.bdiagPreco.L(:,:,k)\...
      [z(istart1:iend1);z(istart2:iend2);z(istart3:iend3)]);

  Pz(istart1:iend1) = zblock(1:2*N);
  Pz(istart2:iend2) = zblock(2*N+1:2*N+2);
  Pz(istart3:iend3) = zblock(2*N+3:2*N+3);
end


end % preconditioner


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vInf = bgFlow(o,X,options)

N = size(X,1)/2;
nv = size(X,2);
[x,y] = oc.getXY(X);

if options.confined
    switch options.farField
        case 'shear'
          vInf = [X(end/2+1:end,:);zeros(N,nv)];  
        
        case 'extenstional'
          vInf = [X(1:end/2,:);-X(end/2+1:end,:)];
          
        otherwise
           vInf = [ones(N,nv);zeros(N,nv)];           
    end
else
    switch options.farField
        case 'constant'
          vInf = [ones(N,nv);zeros(N,nv)];    
          
        case 'couette'
            vInf = [zeros(2*N,1) 1*[-y(:,2)+mean(y(:,2));x(:,2)-mean(x(:,2))]];
            
        otherwise
         vInf = [zeros(N,nv);zeros(N,nv)];     
    end
end

end % bgFlow


end % methods

end % classdef


