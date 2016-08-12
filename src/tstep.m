%TODO: Just realized that o.D and preconditioner do not changes since
%the geometry only rotates and translates.  So, they only need to be
%computed once at the beginning of the code and used at all time steps!
classdef tstep < handle
% This class defines the functions required to advance the geometry
% forward in time.  Routines that we may need to add are different
% integrators such as semi-implicit, SDC, Runge-Kutta, adaptivity


properties
XCoeff     % coefficients for forming the configuration where the
           % layer potential operators are constructed.  Comes from IMEX
rhsCoeff   % coefficients for forming the right-hand side used in IMEX
beta       % single real number used in IMEX
order      % time stepping order (will just do IMEX-Euler for now)
dt         % time step size
DFibers    % Stokes double-layer potential for fibers
DWalls     % Stokes double-layer potential for walls
DupFibers  % Upsampled Stokes double-layer potential matrix for fibers
DupWalls   % Upsampled Stokes double-layer potential matrix for walls
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(options, prams, om, geom, walls)
% o.tstep(options,prams): constructor.  Initialize class.  Take all
% elements of options and prams needed by the time stepper

N = geom.N;
nv = geom.nv;

o.order = options.tstep_order;  % times stepping order (1 for now)
o.dt = prams.T/prams.m;   % time step size

o.ifmm = options.ifmm;   % fast multipole method
o.inear = options.inear; % near-singular integration interpolation scheme
o.gmresTol = prams.gmresTol;
o.farField = @(X) o.bgFlow(X,options.farField);
o.usePreco = options.usePreco;
o.opFibers = poten(prams.N, om);

if options.confined
    o.opWall = poten(prams.Nbd, om);
else
    o.opWall = [];
end

o.profile = options.profile;
o.om = om;

o.DFibers = o.opFibers.stokesDLmatrix(geom);
o.DWalls = o.opWall.stokesDLmatrix(walls);

% create upsampled matrices
Xsou = geom.X; % source positions
Nsou = size(Xsou,1)/2; % number of source points
Nup = Nsou*ceil(sqrt(Nsou));

Xup = [interpft(Xsou(1:Nsou,:),Nup);...
       interpft(Xsou(Nsou+1:2*Nsou,:),Nup)];

geomUp = capsules([],Xup);

o.DupFibers = o.opFibers.stokesDLmatrix(geomUp);

Xsou = walls.X; % source positions
Nsou = size(Xsou,1)/2; % number of source points
Nup = Nsou*ceil(sqrt(Nsou));

Xup = [interpft(Xsou(1:Nsou,:),Nup);...
       interpft(Xsou(Nsou+1:2*Nsou,:),Nup)];

wallsUp = capsules([],Xup);

o.DupWalls = o.opFibers.stokesDLmatrix(wallsUp);


if o.usePreco
    
    if o.profile
        tic;
    end
    
  o.bdiagPreco.L = zeros(2*N+3,2*N+3,nv);
  o.bdiagPreco.U = zeros(2*N+3,2*N+3,nv);
  for k = 1:nv
    [o.bdiagPreco.L(:,:,k),o.bdiagPreco.U(:,:,k)] = lu([...
      -1/2*eye(2*N)+o.D(:,:,k) ...
        [-ones(N,1);zeros(N,1)] ...
        [zeros(N,1);-ones(N,1)] ...
        [geom.X(end/2+1:end,k)-geom.center(2,k);...
          -geom.X(1:end/2,k)+geom.center(1,k)];
        [-geom.sa(:,k)'*2*pi/N zeros(1,N) 0 0 0];
        [zeros(1,N) -geom.sa(:,k)'*2*pi/N 0 0 0];
        [(-geom.X(end/2+1:end,k)'+geom.center(2,k)).*...
            geom.sa(:,k)'*2*pi/N ...
         (geom.X(1:end/2,k)'-geom.center(1,k)).*...
            geom.sa(:,k)'*2*pi/N 0 0 0]]);
    % whole diagonal block which doesn't have a null space
  end
  
  if o.profile
      o.om.writeMessage(['Building preconditioner ... ', num2str(toc)]);
  end
end
  
end % constructor: tstep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eta,Up,wp,iter,iflag,res] = timeStep(o,geom,tau)
% [X,iter,iflag] = timeStep(Xstore) takes the current configuration
% in Xstore (can be a three-dimensional array corresponding to  previous
% time steps if doing multistep) and returns a new shape X, the number
% of GMRES iterations, and the GMRES flag 

N = geom.N;
nv = geom.nv;

Nup = N*ceil(sqrt(N));

if o.inear
    if o.profile
        tic;
    end
    
    o.nearStruct = geom.getZone([],1);
    
    if o.profile
        o.om.writeMessage(['getZone ... ', num2str(toc)]);
    end
end

% rotate DLP matrices and preconditioner matrices
for i = 1:nv
    R = spdiags([sin(tau(i))*ones(2*N,1), cos(tau(i))*ones(2*N,1)...
                -sin(tau(i))*ones(2*N,1)], [-N, 0, N], zeros(2*N, 2*N));
    
    Rup = spdiags([sin(tau(i))*ones(2*Nup,1), cos(tau(i))*ones(2*Nup,1)...
                -sin(tau(i))*ones(2*Nup,1)], [-Nup, 0, Nup], zeros(2*Nup, 2*Nup));  
            
    o.D(:,:,i) = R*o.D(:,:,i)*R';
    o.bdiagPreco.L(1:2*N,1:2*N,i) = R*o.bdiagPreco.L(1:2*N,1:2*N,i)*R';
    o.bdiagPreco.U(1:2*N,1:2*N,i) = R*o.bdiagPreco.U(1:2*N,1:2*N,i)*R';
    
    o.Dup(:,:,i) = Rup*o.Dup(:,:,i)*Rup';
end

%op = o.op;

% build double-layer potential matrix for each rigid body
% o.D = opFibers.stokesDLmatrix(geom);

% compute preconditioner if needed
% if o.usePreco
%     
%     if o.profile
%         tic;
%     end
%     
%   o.bdiagPreco.L = zeros(2*N+3,2*N+3,nv);
%   o.bdiagPreco.U = zeros(2*N+3,2*N+3,nv);
%   for k = 1:nv
%     [o.bdiagPreco.L(:,:,k),o.bdiagPreco.U(:,:,k)] = lu([...
%       -1/2*eye(2*N)+o.D(:,:,k) ...
%         [-ones(N,1);zeros(N,1)] ...
%         [zeros(N,1);-ones(N,1)] ...
%         [geom.X(end/2+1:end,k)-geom.center(2,k);...
%           -geom.X(1:end/2,k)+geom.center(1,k)];
%         [-geom.sa(:,k)'*2*pi/N zeros(1,N) 0 0 0];
%         [zeros(1,N) -geom.sa(:,k)'*2*pi/N 0 0 0];
%         [(-geom.X(end/2+1:end,k)'+geom.center(2,k)).*...
%             geom.sa(:,k)'*2*pi/N ...
%          (geom.X(1:end/2,k)'-geom.center(1,k)).*...
%             geom.sa(:,k)'*2*pi/N 0 0 0]]);
%     % whole diagonal block which doesn't have a null space
%   end
%   
%   if o.profile
%       o.om.writeMessage(['Building preconditioner ... ', num2str(toc)]);
%   end
% end

rhs = o.farField(geom.X);
rhs = -rhs(:);
rhs = [rhs; zeros(3*nv,1)];
% right-hand side with an arbitrarily chosen background flow
% need to negate right-hand side since it moves to the other side of the
% governing equations

maxit = 2*N*nv;%should be a lot lower than this
if ~o.usePreco
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom),rhs,[],o.gmresTol,maxit);
else
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom),rhs,[],o.gmresTol,maxit,@o.preconditionerBD);
end
% Use GMRES to find density function, translation, and rotation
% velocities
iter = I(2);

eta = zeros(2*N,nv);
for k = 1:nv
  eta(1:2*N,k) = Xn((k-1)*2*N+1:k*2*N);
end
% take column vector coming out of GMRES and organize in matrix where
% each column corresponds to the density function of the rigid body

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

N = geom.N; % points per body
nv = geom.nv; % number of bodies

potFibers = o.opFibers;
potWalls = o.opWall;

if o.confined
    Nbd = walls.N;
    nvbd = walls.nv;
else
    Nbd = 0;
    nv = 0;
end

% Output of Tx that corresponds to the velocity of the fibers
valFibers = zeros(2*N,nv);
% Output of Tx that corresponds to the velocity of the walls
valWalls = zeros(2*Nbd,nvbd);
% output of Tx that corresponds to force caused by density function
valForce = zeros(2,nv);
% output of Tx that corresponds to torque caused by density function
valTorque = zeros(nv,1);

% extract density function for fibers
etaFibers = zeros(2*N,nv);
for k = 1:nv
  etaFibers(:,k) = Xn((k-1)*2*N+1:k*2*N);
end

% extract density function for walls
etaWalls = zeros(2*Nbd,nvbd);
for k = 1:nvbd
  etaWalls(:,k) = Xn(2*N*nv+1+(k-1)*Nbd:2*N*nv+1+k*Nbd);
end

% extract translational and rotational velocities
Up = zeros(2,nv);
for k = 1:nv
  Up(:,k) = Xn(2*N*nv+2*Nbd*nvbd+1+2*(k-1):2*N*nv+2*Nbd*nvbd+1+2*k);
end
wp = zeros(1,nv);
for k = 1:nv
  wp(k) = Xn(2*N*nv+2*Nbd*nvbd+2*nv+k);
end

% Jump term in double-layer potential
valFibers = valFibers - 1/2*etaFibers;
valWalls = valWalls - 1/2*etaWalls;

% self contribution
valFibers = valFibers + potFibers.exactStokesDLdiag(geom, o.DFibers, etaFibers);
valWalls = valWalls + potWalls.exactStokesDLdiag(walls, o.DWalls, etaWalls);


% START COMPUTING FIBRE-FIBRE DLP

% kernel function for evaluating the double layer potential.  kernel
% can be a call to a FMM, but for certain computations, direct
% summation is faster, so also want kernelDirect
if o.ifmm
  kernel = @potFibers.exactStokesDLfmm;
else
  kernel = @potFibers.exactStokesDL;
end

kernelDirect = @potFibers.exactStokesDL;

if o.inear
  DLP = @(X) potFibers.exactStokesDLdiag(geom,o.DFibers,X) - 1/2*X;
  ffdlp = potFibers.nearSingInt(geom,eta,DLP,o.nearStruct,kernel,...
        kernelDirect,geom,true,false);
else
  ffdlp = kernel(geom,etaFibers);
end
% END COMPUTING FIBRE-FIBRE DLP


% START COMPUTING WALL-FIBRE DLP
if o.confined
   
   if o.ifmm
        kernel = @potWall.exactStokesDLfmm;       
   else
       kernel = @potwall.exactStokesDL;
   end
   
   kernelDirect = @potWall.exactStokesDL;
   
   if o.inear
       DLP = @(X) potWalls.exactStokesDLdiag(walls, o.DWalls, X) - 1/2*X;
       wfdlp = potWalls.nearSingInt(geom, etaWalls, DLP, o.nearStruct, ...
           kernel, kernelDirect,walls,false,false);
   else
       wfdlp = kernal(geom, etaWalls);
   end
else
    wfdlp = zeros(2*N,nv);
end
% END COMPUTING WALL-FIBRE DLP

% START COMPUTING FIBRE-WALL DLP
if o.confined
   
   if o.ifmm
        kernel = @potWall.exactStokesDLfmm;       
   else
       kernel = @potwall.exactStokesDL;
   end
   
   kernelDirect = @potWall.exactStokesDL;
   
   if o.inear
       DLP = @(X) potWalls.exactStokesDLdiag(walls, o.DWalls, X) - 1/2*X;
       fwdlp = potWalls.nearSingInt(walls, etaWalls, DLP, o.nearStruct, ...
           kernel, kernelDirect,geom,false,false);
   else
       fwdlp = kernal(walls, etaFibers);
   end
else
    fwdlp = zeros(2*N,nv);
end
% END COMPUTING FIBRE-WALL DLP

% EVALUATE VELOCITY ON WALLS
if o.confined
    valWalls = valWalls - 1/2*etaWalls + potWalls.exactStokesDLdiag(walls, o.DWall, etaWalls);
    valWalls(:,1) = valWalls(:,1) + potWalls.exactStokesN0diag(walls, o.N0wall, etaWalls(:,1));
    
    valWalls = valWalls + fwdlp;
end

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

for k = 1:nv
  valForce(1,k) = sum(eta(1:N,k).*geom.sa(:,k))*2*pi/N;
  valForce(2,k) = sum(eta(N+1:2*N,k).*geom.sa(:,k))*2*pi/N;
end

for k = 1:nv
  valTorque(k) = sum(((geom.X(N+1:2*N,k)-geom.center(2,k)).*eta(1:N,k) - ...
                     ((geom.X(1:N,k)-geom.center(1,k)).*eta(N+1:2*N,k))).*...
                     geom.sa(:,k))*2*pi/N;
end

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
function vInf = bgFlow(o,X,type)

N = size(X,1)/2;
nv = size(X,2);

if strcmp(type,'shear')
  vInf = [X(end/2+1:end,:);zeros(N,nv)];
elseif strcmp(type,'extensional')
  vInf = [X(1:end/2,:);-X(end/2+1:end,:)];
elseif strcmp(type, 'poiseuille')
    %poiseuille flow with 0s at -100 and 100
    vInf = [-(X(end/2+1:end,:) - 15).*(X(end/2+1:end,:) + 15)/(10^2); zeros(N,nv)];
elseif strcmp(type, 'rotlet')
    %rotlet centered at (0,0)
    rot = @(x, y, xc, yc) 1./((x-xc).^2+(y-yc).^2).*[(y-yc), -(x-xc)];
    
    vInf = zeros(size(X));
    
    for i = 1:N
        
        for j = 1:nv
            rotTmp = rot(X(i,j), X(i+N,j), 0, 0) + rot(X(i,j), X(i+N,j), 4, 4);
            vInf(i,j) = rotTmp(1);
            vInf(i+N,j) = rotTmp(2);
        end
    end
else
  vInf = [ones(N,nv);zeros(N,nv)];
end

end % bgFlow


end % methods

end % classdef


