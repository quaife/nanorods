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
D          % Stokes double-layer potential matrix
ifmm       % flag for using the FMM
inear      % flag for using near-singular integration
gmresTol   % GMRES tolerance
nearStruct % near-singular integration structure 
farField   % background flow
usePreco   % use a block-diagonal preconditioner
bdiagPreco % block-diagonal preconditioner
op         % class for layer potentials
end % properties


methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(options,prams)
% o.tstep(options,prams): constructor.  Initialize class.  Take all
% elements of options and prams needed by the time stepper

o.order = options.order;  % times stepping order (1 for now)
o.dt = prams.T/prams.m;   % time step size
[o.XCoeff,o.rhsCoeff,o.beta] = o.getCoeff(o.order);
% load variables for IMEX integrator
o.ifmm = options.ifmm;
o.inear = options.inear; % near-singular integration interpolation scheme
o.gmresTol = prams.gmresTol;
o.farField = @(X) o.bgFlow(X,options.farField);
o.usePreco = options.usePreco;
o.op = poten(prams.N);

end % constructor: tstep



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eta,Up,wp,iter,iflag,res] = timeStep(o,geom)
% [X,iter,iflag] = timeStep(Xstore) takes the current configuration
% in Xstore (can be a three-dimensional array corresponding to  previous
% time steps if doing multistep) and returns a new shape X, the number
% of GMRES iterations, and the GMRES flag 

N = geom.N;
nv = geom.nv;

Xm = geom.X;
% this'll be a linear comination of Xstore if doing multistep

if o.inear
    o.nearStruct = geom.getZone([],1);
end

op = o.op;

% build double-layer potential matrix for each rigid body
o.D = op.stokesDLmatrix(geom);

% compute preconditioner if needed
if o.usePreco
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
end

%Z1 = zeros(nv*(2*N+3));
%Z2 = zeros(nv*(2*N+3));
%for k = 1:nv*(2*N+3);
%  e = zeros(nv*(2*N+3),1); e(k) = 1;
%  Z1(:,k) = o.timeMatVec(e,geom);
%  Z2(:,k) = o.preconditionerBD(e);
%end

%figure(1);
%spy(Z1)
%figure(2);
%spy(inv(Z2))
%figure(3);
%Z2Z1 = Z2*Z1;
%surf(Z2Z1 - eye(2*N*nv+3*nv))
%pause

rhs = [o.farField(geom.X)];
%rhs = [geom.X(end/2+1:end,:);zeros(N,nv);zeros(3,nv)];
rhs = -rhs(:);
rhs = [rhs; zeros(3*nv,1)];
% right-hand side with an arbitrarily chosen background flow
% need to negate right-hand side since it moves to the other side of the
% governing equations

%%z = o.timeMatVec(rhs,geom);
%z = o.preconditionerBD(rhs);
%for k = 1:2*nv
%  is = (k-1)*N+1;
%  ie = is + N - 1;
%%  semilogy(abs(fftshift(fft(z(is:ie))))/N)
%  plot(z(is:ie))
%  pause
%end
%%plot(o.preconditionerBD(o.timeMatVec(rhs,geom)) -rhs);
%%pause

maxit = 2*N*nv;%should be a lot lower than this
if ~o.usePreco
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom),...
      rhs,[],o.gmresTol,maxit);
else
  [Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom),...
      rhs,[],o.gmresTol,maxit,@o.preconditionerBD);
end
% Use GMRES to find density function, translation, and rotation
% velocities
iter = I(2);

oc = curve;
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
function Tx = timeMatVec(o,Xn,geom)
% Tx = timeMatVec(Xn,geom) does a matvec for GMRES with the IMEX scheme
% Xn is a column vector which contains the x-coordinate followed by the
% y-coordinate of body 1, x-coordinate followed by the y-coordinate of
% body 2, etc

N = geom.N; % points per body
op = o.op; 
% this can be setup as a property of tstep if it is too expensive
nv = geom.nv; % number of bodies

valPos = zeros(2*N,nv);
% Output of Tx that corresponds to the position of the bodies
valForce = zeros(2,nv);
% output of Tx that corresponds to force caused by density function
valTorque = zeros(nv,1);
% output of Tx that corresponds to torque caused by density function

eta = zeros(2*N,nv);
for k = 1:nv
  eta(1:2*N,k) = Xn((k-1)*2*N+1:k*2*N);
  %disp((k-1)*2*N+1:k*2*N)
end
% organize as matrix where columns correspond to unique bodies
Up = zeros(2,nv);
for k = 1:nv
  Up(:,k) = Xn(2*N*nv+(k-1)*2+1:2*N*nv+k*2);
end
wp = zeros(1,nv);
for k = 1:nv
  wp(k) = Xn(2*N*nv+2*nv+k);
end

valPos = valPos - 1/2*eta;
% Jump term in double-layer potential

valPos = valPos + op.exactStokesDLdiag(geom,o.D,eta);
% self contribution


%theta = (0:63)'*2*pi/64;
%for k = 1:49
%  eta(:,k) = [exp(cos(theta));cos(theta)];
%end
%eta = [[exp(cos(theta));cos(theta)] [sin(cos(theta));ones(128,1)]];

if o.ifmm
  kernel = @op.exactStokesDLfmm;
else
  kernel = @op.exactStokesDL;
end
kernelDirect = @op.exactStokesDL;
% kernel function for evaluating the double layer potential.  kernel
% can be a call to a FMM, but for certain computations, direct
% summation is faster, so also want kernelDirect
if o.inear
  DLP = @(X) op.exactStokesDLdiag(geom,o.D,X) - 1/2*X;
  Fdlp = op.nearSingInt(geom,eta,DLP,...
      o.nearStruct,kernel,kernelDirect,geom,true,false);
else
  Fdlp = kernel(geom,eta);
end
valPos = valPos + Fdlp;
% Add in contribution from other bodies using exactStokesDL

for k = 1:nv
  valPos(1:end/2,k) = valPos(1:end/2,k) - Up(1,k);
  valPos(end/2+1:end,k) = valPos(end/2+1:end,k) - Up(2,k);
end

for k = 1:nv
  valPos(1:end/2,k) = valPos(1:end/2,k) + ...
      (geom.X(end/2+1:end,k) - geom.center(2,k))*wp(k);
  valPos(end/2+1:end,k) = valPos(end/2+1:end,k) - ...
      (geom.X(1:end/2,k) - geom.center(1,k))*wp(k);
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

Tx = [valPos(:);-valForce(:);-valTorque(:)];

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
  vInf = [-X(1:end/2,:);X(end/2+1:end,:)];
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XCoeff,rhsCoeff,beta] = getCoeff(o,order)
% [XCoeff,rhsCoeff,beta] = getCoeff(order) generates the coefficients
% required to discretize the derivative.  First-order time  derivatives
% are approximated by beta*x^{N+1} + rhscoeff.*[x^{N} x^{N-1} ...]
% Explicit terms (operators) are discretized at XCoeff.*[x^{N} x^{N-1}
% ...] All rules are from Ascher, Ruuth, Wetton 1995.

beta = 1;
XCoeff = 1;
rhsCoeff = 1;


end % getCoeff


end % methods

end % classdef


