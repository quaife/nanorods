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
o.confined = options.confined;

o.dt = prams.T/prams.m;                  
o.gmresTol = prams.gmresTol;             

o.farField = @(X) o.bgFlow(X,options); 
o.om = om;

N = prams.N;
Nbd = prams.Nbd;
nv = prams.nv;
nbd = prams.nbd;

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

% CREATE UPSAMPLED MATRICES
% FIBRE-FIBRE
Xsou = geom.X; 
Nup = N*ceil(sqrt(N));

Xup = [interpft(Xsou(1:N,:),Nup);...
       interpft(Xsou(N+1:2*N,:),Nup)];

geomUp = capsules([],Xup);
o.Dupf = o.opFibers.stokesDLmatrix(geomUp);

% WALL-WALL
Xsou = walls.X; 
Nup = Nbd*ceil(sqrt(Nbd));

Xup = [interpft(Xsou(1:Nbd,:),Nup);...
       interpft(Xsou(Nbd+1:2*Nbd,:),Nup)];

wallsUp = capsules([],Xup);
o.Dupw = o.opFibers.stokesDLmatrix(wallsUp);

% CREATE N0 MATRIX
o.N0w = o.opWall.stokesN0matrix(walls);

% CREATE BLOCK-DIAGONAL PRECONDITIONER
if o.usePreco
    
    if o.profile
        tic;
    end
    
    if o.confined
        o.precoF.L = zeros(2*N+3 + 2*Nbd,2*N+3 + 2*Nbd,nv);
        o.precoF.U = zeros(2*N+3 + 2*Nbd,2*N+3 + 2*Nbd,nv);
        for k = 1:nv
            [o.precoF.L(:,:,k),o.precoF.U(:,:,k)] =...
                lu([-1/2*eye(2*N)+o.Df(:,:,k) ...
                zeros(2*N, 2*Nbd)...
                [-ones(N,1);zeros(N,1)] ...
                [zeros(N,1);-ones(N,1)] ...
                [geom.X(end/2+1:end,k)-geom.center(2,k); -geom.X(1:end/2,k)+geom.center(1,k)];...
                [-geom.sa(:,k)'*2*pi/N zeros(1,N + 2*Nbd) 0 0 0];
                [zeros(1,N + 2*Nbd) -geom.sa(:,k)'*2*pi/N 0 0 0];
                [(-geom.X(end/2+1:end,k)'+geom.center(2,k)).*geom.sa(:,k)'*2*pi/N ...
                (geom.X(1:end/2,k)'-geom.center(1,k)).*geom.sa(:,k)'*2*pi/N zeros(1, 2*Nbd + 3)]]);
        end
        
        o.precoF.L = zeros(2*N+3 + 2*Nbd,2*N+3 + 2*Nbd,nv);
        o.precoF.U = zeros(2*N+3 + 2*Nbd,2*N+3 + 2*Nbd,nv);
    else
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
    o.rhs = [-rhs; zeros(2*Nbd*nbd,1); zeros(3*nv,1)];
end

end % constructor: tstep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [etaF,etaW,Up,wp,iter,iflag,res] = timeStep(o,geom,tau,walls)
% [X,iter,iflag] = timeStep(Xstore) takes the current configuration
% in Xstore (can be a three-dimensional array corresponding to  previous
% time steps if doing multistep) and returns a new shape X, the number
% of GMRES iterations, and the GMRES flag 

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

Nbd = walls.N;
nbd = walls.nv;

Nup = N*ceil(sqrt(N));

for i = 1:nv
    R = spdiags([sin(tau(i))*ones(2*N,1), cos(tau(i))*ones(2*N,1)...
                -sin(tau(i))*ones(2*N,1)], [-N, 0, N], zeros(2*N, 2*N));
    
    Rup = spdiags([sin(tau(i))*ones(2*Nup,1), cos(tau(i))*ones(2*Nup,1)...
                -sin(tau(i))*ones(2*Nup,1)], [-Nup, 0, Nup], zeros(2*Nup, 2*Nup));  
            
    o.Df(:,:,i) = R*o.Df(:,:,i)*R';
    o.Dupf(:,:,i) = Rup*o.Dupf(:,:,i)*Rup';
    
    if o.usePreco
        o.precoF.L(1:2*N,1:2*N,i) = R*o.bdiagPreco.L(1:2*N,1:2*N,i)*R';
        o.precoF.U(1:2*N,1:2*N,i) = R*o.bdiagPreco.U(1:2*N,1:2*N,i)*R;
    end    
end

% SOLVE SYSTEM USING GMRES
maxit = 2*N*nv;% should be a lot lower than this
if ~o.usePreco
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom,walls),o.rhs,[],o.gmresTol,...
      maxit);
else
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom,walls),o.rhs,[],o.gmresTol,...
      maxit,@o.preconditionerBD);
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
etaW = zeros(2*Nbd, nbd);
for k = 1:nbd
   etaW(:,k) = Xn(2*N+1+(k-1)*2*Nbd:2*N+k*2*Nbd);
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

end % timeStep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tx = timeMatVec(o,Xn,geom,walls)
% Tx = timeMatVec(Xn,geom) does a matvec for GMRES 
% Xn is a state vector that contains the following in order:
% fiber1: eta_x, eta_y fiber2: eta_x, eta_y ... fiberNv: eta_x, eta_y
% wall1 : xi_x, xi_y wall2: xi_x, xi_y ... wallnbd: xi_x, xi_y
% fiber1: u, v fiber2: u, v ... fiberNv: u, v
% fiber1: omega fiber2: omega ... fiberNv: omega

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
    nv = 0;
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
for k = 1:nbd
  etaW(:,k) = Xn(2*N*nv+1+(k-1)*2*Nbd:2*N*nv+k*2*Nbd);
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
    xi = zeros(2,nbd-1);
    for k = 1:nbd-1
       xi(:,k) = Xn(2*N*nv+2*Nbd*nbd+3*nv+1+2*(k-1):2*N*nv+2*Nbd*nbd+3*nv+2*k); 
    end

    lambda=zeros(1,nbd-1);
    for k=1:nbd-1
        lambda(k) = Xn(2*N*nv+2*Nbd*nbd+3*nv+2*(nbd-1)+k);
    end    
else
    xi = zeros(2,1);
    lambda = 0;
end

% END FORMATTING UNKNOWN VECTOR

% ADD JUMP IN DLP
valFibers = valFibers - 1/2*etaF;
valWalls = valWalls - 1/2*etaW;

% ADD SELF CONTRIBUTION
valFibers = valFibers + potFibers.exactStokesDLdiag(geom, o.Df, etaF);
valWalls = valWalls + potWalls.exactStokesDLdiag(walls, o.Dw, etaW);

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
  ffdlp = kernel(geom, etaFibers);
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
    wfdlp = zeros(2*N,nbd);
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
    fwdlp = zeros(2*Nbd,nv);
end
% END OF TARGET == FIBRES

% START OF TARGET == WALLS
if o.confined
    
    if o.ifmm
        kernel = @potWalls.exactStokesDLfmm;
    else
        kernel = @potwalls.exactStokesDL;
    end
    
    wwdlp = kernel(walls, etaW);
    
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
        fsr = fsr + o.RSlets(geom.X, walls.center(:,k), xi(:,k-1), lambda(k-1));
    end 
    
    % END TARGETS == FIBRES
    
    % START TARGETS == WALLS
    wsr = 0;
    for k = 2:nbd % loop over all walls, except outer wall
        wsr = wsr + o.RSlets(walls.X, walls.center(:,k), xi(:,k-1), lambda(k-1));
    end   
    
    % END TARGETS == WALLS
    % END SOURCE == ROTLETS
    
    % EVALUATE ROTLET AND STOKESLET EQUATIONS
    z = o.letsIntegrals([xi;lambda], etaW, walls);
    valStokeslets = z(1:2*(nbd-1));
    valRotlets = z(2*(nbd-1)+1:end);
else
    fsr = 0;
    wsr = 0;
    valRotlets = [];
    valStokeslets = [];
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
Tx = [valFibers(:); valWalls(:); -valForce(:);-valTorque(:);valRotlets(:);valStokeslets(:)];

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

        case 'couette'
            %[-y,x]==>clockwise
            vInf = 1*[-y(:,1)+mean(y(:,1));x(:,1)-mean(x(:,1))];
            %[y,-x]==>counter clockwise
            vInf = [vInf; 4*[y(:,2)+mean(y(:,2)); -x(:,2)-mean(x(:,2))]];

        otherwise
            vInf = [zeros(N,nv);zeros(N,nv)];
    end
else
    switch options.farField
        case 'shear'
            vInf = [X(end/2+1:end,:);zeros(N,nv)];

        case 'extenstional'
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
       M(:,k) = o.timeMatVec(I(:,k),geom,walls) + 0.5*I(:,k);
    end    
end % build_matrix
end % methods

end % classdef


