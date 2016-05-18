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
o.ifmm = false; % TODO: will need this to be true later
o.inear = options.inear; % near-singular integration interpolation scheme
o.gmresTol = prams.gmresTol;
o.farField = @(X) o.bgFlow(X,options.farField);
o.op = poten(prams.N);


end % constructor: tstep



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eta,Up,wp,iter,iflag] = timeStep(o,geom)
% [X,iter,iflag] = timeStep(Xstore) takes the current configuration
% in Xstore (can be a three-dimensional array corresponding to  previous
% time steps if doing multistep) and returns a new shape X, the number
% of GMRES iterations, and the GMRES flag 

N = geom.N;
nv = geom.nv;

Xm = geom.X;
% this'll be a linear comination of Xstore if doing multistep

o.nearStruct = geom.getZone([],1);

op = o.op;
o.D = op.stokesDLmatrix(geom);
% build double-layer potential matrix for each rigid body

rhs = [o.farField(geom.X);zeros(3,nv)];
%rhs = [geom.X(end/2+1:end,:);zeros(N,nv);zeros(3,nv)];
rhs = -rhs(:);
% right-hand side with an arbitrarily chosen background flow
% need to negate right-hand side since it moves to the other side of the
% governing equations

maxit = 2*N;
[Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom),...
    rhs,[],o.gmresTol,maxit);
% Use GMRES to find new geometry
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


kernel = @op.exactStokesDL;
kernelDirect = @op.exactStokesDL;
% kernel function for evaluating the double layer potential.  kernel
% can be a call to a FMM, but for certain computations, direct
% summation is faster, so also want kernelDirect
%theta = (0:127)'*2*pi/128;
%eta = [[exp(cos(theta));cos(theta)] [sin(cos(theta));ones(128,1)]];
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
function vInf = bgFlow(o,X,type)

N = size(X,1)/2;
nv = size(X,2);

if strcmp(type,'shear')
  vInf = [X(end/2+1:end,:);zeros(N,nv)];
elseif strcmp(type,'extensional')
  vInf = [-X(1:end/2,:);X(end/2+1:end,:)];
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


