classdef poten 
% this class defines single and double layers for various kernels
% (stokes, laplace) on 2D periodic curves.
% Defines the matricies that map a density function defined on the
% boundary of a curve to the layer potential evaluated on the curve,
% and defines the operators that takes a density function defined on
% the boundary of a curve and returns the layer 
% potential at arbitrary target points.
% This class also has the main routine that evaluates layer
% potentials using near-singular integration.
    
properties
  N; % points per curve
  interpMat;  
  % upsampled quadrature rules for Alpert's quadrature rule.
  gamma
  % indicies for Kapur-Rokhlin quadrature

end % properties

methods 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = poten(N)
% o = poten(N): constructor; N is the number of points per curve

o.N = N;
o.interpMat = o.lagrangeInterp;
% load in the interpolation matrix which is precomputed with 7
% interpolation points

%o.gamma = [ +1.825748064736159e0;...
%            -1.325748064736161e0];
% second-order accurate
o.gamma = [+4.967362978287758e+0; ...
         -1.620501504859126e+1; ...
         +2.585153761832639e+1; ...
         -2.222599466791883e+1; ...
         +9.930104998037539e+0; ...
         -1.817995878141594e+0];
% sixth-order accurate
% Load weights that are required near the log singularity

end % poten: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LP = nearSingInt(o,vesicleSou,f,selfMat,...
    NearStruct,kernel,kernelDirect,vesicleTar,tEqualS,upNear,idebug)
% LP = nearSingInt(vesicle,f,selfMat,NearStruct,kernel,kernelDirect,
% vesicleTar,tEqualS,idebug) computes a layer potential due to f at all
% points in vesicleTar.X.  If tEqualS==true, then the vesicleTar ==
% vesicleSou and the self-vesicle interaction is skipped.  selfMat is
% the diagonal of the potential needed to compute the layer potential of
% each vesicle indepenedent of all others.  kernel and kernelDirect are
% two (possibly the same) routines that compute the layer potential.
% kernelDirect always uses the direct method whereas kernel may use an
% FMM-accelerated method.  NearStruct is a structure containing the
% variables zone,dist,nearest,icp,argnear which are required by
% near-singular integration (they keep everything sorted and
% precomputed) Everything is in the 2*N x nv format Can pass a final
% argument if desired so that plots of the near-singular integration
% algorithm are displayed

if (tEqualS && size(vesicleSou.X,2) == 1)
  LP = zeros(size(vesicleSou.X));
  return
end
% only a single vesicle, so velocity on all other vesicles will always
% be zero

if (nargin == 10)
  idebug = false;
end

dist = NearStruct.dist;
zone = NearStruct.zone;
nearest = NearStruct.nearest;
icp = NearStruct.icp;
argnear = NearStruct.argnear;

Xsou = vesicleSou.X; % source positions
Nsou = size(Xsou,1)/2; % number of source points
nvSou = size(Xsou,2); % number of source 'vesicles'
Xtar = vesicleTar.X; % target positions
Ntar = size(Xtar,1)/2; % number of target points
nvTar = size(Xtar,2); % number of target 'vesicles'

h = vesicleSou.length/Nsou; % arclength term

Nup = Nsou*ceil(sqrt(Nsou));

vself = selfMat(f);

% upsample to N^(3/2).  
Xup = [interpft(Xsou(1:Nsou,:),Nup);...
       interpft(Xsou(Nsou+1:2*Nsou,:),Nup)];
fup = [interpft(f(1:Nsou,:),Nup);...
       interpft(f(Nsou+1:2*Nsou,:),Nup)];

vesicleUp = capsules(Xup,[],[],...
    vesicleSou.kappa,vesicleSou.viscCont,vesicleSou.antiAlias);
% Build an object with the upsampled vesicle

interpOrder = size(o.interpMat,1);
% lagrange interpolation order
p = ceil((interpOrder+1)/2);
% want half of the lagrange interpolation points to the left of the
% closest point and the other half to the right

if tEqualS % sources == targets
  if nvSou > 1
    if (strfind(char(kernel),'fmm'))
      farField = kernel(vesicleUp,fup);
      farField = farField(1:Nup/Ntar:end,:);
      % evaluate layer potential at all targets except ignore the
      % diagonal term
    else
      for k = 1:nvSou
        K = [(1:k-1) (k+1:nvSou)];
        [~,farField(:,k)] = kernelDirect(vesicleUp,fup,Xtar(:,k),K);
      end
      % This is a huge savings if we are using a direct method rather
      % than the fmm to evaluate the layer potential.  The speedup is
      % more than N^{1/2}, where N is the resolution of the vesicles
      % that we are computing with
    end
  else
    farField = zeros(2*Ntar,nvTar);
  end

else % sources ~= targets
  [~,farField] = kernel(vesicleUp,fup,Xtar,1:nvSou);
  % evaluate layer potential due to all 'vesicles' at all points in
  % Xtar
end
% Use upsampled trapezoid rule to compute layer potential

nearField = zeros(2*Ntar,nvTar);

beta = 1.1;
% small buffer to make sure Lagrange interpolation points are not in the
% near zone
for k1 = 1:nvSou
  if tEqualS % sources == targets
    K = [(1:k1-1) (k1+1:nvTar)];
    % skip diagonal vesicle
  else % sources ~= targets
    K = (1:nvTar);
    % consider all vesicles
  end
  for k2 = K
    J = find(zone{k1}(:,k2) == 1);
    % set of points on vesicle k2 close to vesicle k1
    if (numel(J) ~= 0)
      indcp = icp{k1}(J,k2);
      % closest point on vesicle k1 to each point on vesicle k2 
      % that is close to vesicle k1
      for j = 1:numel(J)
          pn = mod((indcp(j)-p+1:indcp(j)-p+interpOrder)' - 1,Nsou) + 1;
          % index of points to the left and right of the closest point
          v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
            o.interpMat*vself(pn,k1));
          vel(J(j),k2,k1) = v(end);  
          % x-component of the velocity at the closest point
          v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
            o.interpMat*vself(pn+Nsou,k1));
          vel(J(j)+Ntar,k2,k1) = v(end);
          % y-component of the velocity at the closest point
      end
%     compute values of velocity at required intermediate points
%     using local interpolant
      
      if ((numel(J) + numel(fup)) >= 512 && numel(J) > 32)
        [~,potTar] = kernel(vesicleUp,fup,...
           [Xtar(J,k2);Xtar(J+Ntar,k2)],k1);
      else
        [~,potTar] = kernelDirect(vesicleUp,fup,...
           [Xtar(J,k2);Xtar(J+Ntar,k2)],k1);
      end
      % Need to subtract off contribution due to vesicle k1 since its
      % layer potential will be evaulted using Lagrange interpolant of
      % nearby points
      nearField(J,k2) =  - potTar(1:numel(J));
      nearField(J+Ntar,k2) =  - potTar(numel(J)+1:end);
      
      XLag = zeros(2*numel(J),interpOrder - 1);
      % initialize space for initial tracer locations
      for i = 1:numel(J)
        nx = (Xtar(J(i),k2) - nearest{k1}(J(i),k2))/...
            dist{k1}(J(i),k2);
        ny = (Xtar(J(i)+Ntar,k2) - nearest{k1}(J(i)+Ntar,k2))/...
            dist{k1}(J(i),k2);
        XLag(i,:) = nearest{k1}(J(i),k2) + ...
            beta*h*nx*(1:interpOrder-1);
        XLag(i+numel(J),:) = nearest{k1}(J(i)+Ntar,k2) + ...
            beta*h*ny*(1:interpOrder-1);
        % Lagrange interpolation points coming off of vesicle k1 All
        % points are behind Xtar(J(i),k2) and are sufficiently far from
        % vesicle k1 so that the Nup-trapezoid rule gives sufficient
        % accuracy
      end

      if (numel(XLag)/2 > 100)
        [~,lagrangePts] = kernel(vesicleUp,fup,XLag,k1);
      else
        [~,lagrangePts] = kernelDirect(vesicleUp,fup,XLag,k1);
      end
      % evaluate velocity at the lagrange interpolation points
      
      for i = 1:numel(J)
        Px = o.interpMat*[vel(J(i),k2,k1) ...
            lagrangePts(i,:)]';
        Py = o.interpMat*[vel(J(i)+Ntar,k2,k1) ...
            lagrangePts(i+numel(J),:)]';
        % Build polynomial interpolant along the one-dimensional
        % points coming out of the vesicle
        dscaled = full(dist{k1}(J(i),k2)/(beta*h*(interpOrder-1)));
        % Point where interpolant needs to be evaluated

        v = filter(1,[1 -dscaled],Px);
        nearField(J(i),k2) = nearField(J(i),k2) + ...
            v(end);
        v = filter(1,[1 -dscaled],Py);
        nearField(J(i)+Ntar,k2) = nearField(J(i)+Ntar,k2) + ...
            v(end);

        if idebug
          figure(2); clf; hold on;
          plot(Xsou(1:Nsou,:),Xsou(Nsou+1:end,:),'r.','markersize',10)
          plot(Xtar(1:Ntar,:),Xtar(Ntar+1:end,:),'k.','markersize',10)
          plot(Xtar(J,k2),Xtar(Ntar+J,k2),'b.','markersize',10)
          plot(XLag(1:numel(J),:),XLag(numel(J)+1:end,:),'kx','markersize',10)
          plot(XLag(i,:),XLag(numel(J)+i,:),'gx','markersize',10)
          axis equal

          figure(1); clf; hold on
          plot((0:interpOrder-1)*beta*h,...
              real([vel(J(i),k2,k1) lagrangePts(i,:)]),'g-o')
          plot((0:interpOrder-1)*beta*h,...
              real([vel(J(i)+Ntar,k2,k1) lagrangePts(i+numel(J),:)]),'r--o')
          pause
        end
        % DEBUG: PASS IN idebug=true INTO THIS ROUTINE AND THEN YOU CAN SEE
        % THE INTERPOLATION POINTS AND CHECK THE SMOOTHNESS OF THE INTERPOLANT

      end % i
    end % numel(J) ~= 0
    % Evaluate layer potential at Lagrange interpolation
    % points if there are any
  end % k2
end % k1
% farField

LP = farField + nearField;
% Add kernel due to far points and near points.  Far points were
% upsampled if source==vesicle so need to truncate here.  We are 
% only using Ntar target points.  Note that it is only the sources 
% that were upsampled

end % nearSingInt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF ROUTINES THAT BUILD LAYER-POTENTIAL MATRICIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = stokesDLmatrix(o,vesicle)
% D = stokesDLmatrix(vesicle), generate double-layer potential for 
% Stokes vesicle is a data structure defined as in the capsules class
% D is (2N,2N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X 

oc = curve;
[x,y] = oc.getXY(vesicle.X);
% Vesicle positions
N = vesicle.N;
% number of points per vesicles
D = zeros(2*N,2*N,vesicle.nv);
% initialize space for double-layer potential matrix

for k=1:vesicle.nv  % Loop over curves
  xx = x(:,k);
  yy = y(:,k);
  % locations
  [tx,ty] = oc.getXY(vesicle.xt);
  tx = tx(:,k); ty = ty(:,k);
  % Vesicle tangent
  sa = vesicle.sa(:,k)';
  % Jacobian
  cur = vesicle.cur(:,k)';
  % curvature

  xtar = xx(:,ones(N,1))';
  ytar = yy(:,ones(N,1))';
  % target points

  xsou = xx(:,ones(N,1));
  ysou = yy(:,ones(N,1));
  % source points

  txsou = tx';
  tysou = ty';
  % tangent at srouces
  sa = sa(ones(N,1),:);
  % Jacobian

  diffx = xtar - xsou;
  diffy = ytar - ysou;
  rho4 = (diffx.^2 + diffy.^2).^(-2);
  rho4(1:N+1:N.^2) = 0;
  % set diagonal terms to 0

  kernel = diffx.*(tysou(ones(N,1),:)) - ...
          diffy.*(txsou(ones(N,1),:));
  kernel = kernel.*rho4.*sa;
  kernel = kernel;

  D11 = kernel.*diffx.^2;
  % (1,1) component
  D11(1:N+1:N.^2) = 0.5*cur.*sa(1,:).*txsou.^2;
  % diagonal limiting term

  D12 = kernel.*diffx.*diffy;
  % (1,2) component
  D12(1:N+1:N.^2) = 0.5*cur.*sa(1,:).*txsou.*tysou;
  % diagonal limiting term

  D22 = kernel.*diffy.^2;
  % (2,2) component
  D22(1:N+1:N.^2) = 0.5*cur.*sa(1,:).*tysou.^2;
  % diagonal limiting term

  D(:,:,k) = [D11 D12; D12 D22];
  % build matrix with four blocks
  D(:,:,k) = 1/pi*D(:,:,k)*2*pi/N;
  % scale with the arclength spacing and divide by pi

end % k

end % stokesDLmatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF ROUTINES THAT BUILD LAYER-POTENTIAL MATRICIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF ROUTINES THAT EVALUATE LAYER-POTENTIALS
% WHEN SOURCES == TARGETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DLP = exactStokesDLdiag(o,vesicle,D,f)
% DLP = exactStokesDLdiag(vesicle,f,K) computes the diagonal term of
% the double-layer potential due to f around all vesicles.  Source and
% target points are the same.  This uses trapezoid rule with the
% curvature at the diagonal in order to guarantee spectral accuracy.
% This routine can either compute the double-layer potential
% matrix-free, which may upsample the number of source points.  Or, if
% the matrix D is passed in and anti-aliasing is not requested, it will
% simply do the matrix-vector product with the precomputed matrix D.

DLP = zeros(2*vesicle.N,vesicle.nv);
for k = 1:vesicle.nv
  DLP(:,k) = D(:,:,k) * f(:,k);
end

end % exactStokesDLdiag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF ROUTINES THAT EVALUATE LAYER-POTENTIALS
% WHEN SOURCES == TARGETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF ROUTINES THAT EVALUATE LAYER-POTENTIALS
% WHEN SOURCES ~= TARGETS.  CAN COMPUTE LAYER POTENTIAL ON EACH
% VESICLE DUE TO ALL OTHER VESICLES (ex. stokesSLP) AND CAN
% COMPUTE LAYER POTENTIAL DUE TO VESICLES INDEXED IN K1 AT 
% TARGET POINTS Xtar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLP,stokesDLPtar] = ...
    exactStokesDL(o,vesicle,f,Xtar,K1)
% [stokesDLP,stokesDLPtar] = exactStokesDL(vesicle,f,Xtar,K1) computes
% the double-layer potential due to f around all vesicles except
% itself.  Also can pass a set of target points Xtar and a collection
% of vesicles K1 and the double-layer potential due to vesicles in K1
% will be evaluated at Xtar.  Everything but Xtar is in the 2*N x nv
% format Xtar is in the 2*Ntar x ncol format
normal = [vesicle.xt(vesicle.N+1:2*vesicle.N,:);...
         -vesicle.xt(1:vesicle.N,:)]; 
% Normal vector

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stokesDLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stokesDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, the user does not need the velocity at arbitrary
  % points
end
den = f.*[vesicle.sa;vesicle.sa]*2*pi/vesicle.N;
% jacobian term and 2*pi/N accounted for here

oc = curve;
[xsou,ysou] = oc.getXY(vesicle.X(:,K1));
xsou = xsou(:); ysou = ysou(:);
xsou = xsou(:,ones(Ntar,1))';
ysou = ysou(:,ones(Ntar,1))';

[denx,deny] = oc.getXY(den(:,K1));
denx = denx(:); deny = deny(:);
denx = denx(:,ones(Ntar,1))';
deny = deny(:,ones(Ntar,1))';

[normalx,normaly] = oc.getXY(normal(:,K1));
normalx = normalx(:); normaly = normaly(:);
normalx = normalx(:,ones(Ntar,1))';
normaly = normaly(:,ones(Ntar,1))';

for k = 1:ncol % loop over columns of target points
  [xtar,ytar] = oc.getXY(Xtar(:,k));
  xtar = xtar(:,ones(vesicle.N*numel(K1),1));
  ytar = ytar(:,ones(vesicle.N*numel(K1),1));
  
  diffx = xtar-xsou; diffy = ytar-ysou;
  dis2 = (diffx).^2 + (diffy).^2;
  % difference of source and target location and distance squared
  
  
  rdotnTIMESrdotf = (diffx.*normalx + diffy.*normaly)./dis2.^2 .* ...
      (diffx.*denx + diffy.*deny);
  % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term
  rdotnTIMESrdotf = rdotnTIMESrdotf * diag(1-vesicle.viscCont(1));
  % have accounted for the scaling with (1-\nu) here
  
  stokesDLPtar(1:Ntar,k) = stokesDLPtar(1:Ntar,k) + ...
      sum(rdotnTIMESrdotf.*diffx,2);
  stokesDLPtar(Ntar+1:end,k) = stokesDLPtar(Ntar+1:end,k) + ...
      sum(rdotnTIMESrdotf.*diffy,2);
  % r \otimes r term of the double-layer potential
end
stokesDLPtar = stokesDLPtar/pi;
% double-layer potential due to vesicles indexed over K1 evaluated at
% arbitrary points

stokesDLP = zeros(2*vesicle.N,vesicle.nv);
if (nargin == 3 && vesicle.nv > 1)
  oc = curve;
  for k = 1:vesicle.nv
    K = [(1:k-1) (k+1:vesicle.nv)];
    [x,y] = oc.getXY(vesicle.X(:,K));
    [nx,ny] = oc.getXY(normal(:,K));
    [denx,deny] = oc.getXY(den(:,K));
    for j=1:vesicle.N
      diffxy = [vesicle.X(j,k) - x ; vesicle.X(j+vesicle.N,k) - y];
      dis2 = diffxy(1:vesicle.N,:).^2 + ...
          diffxy(vesicle.N+1:2*vesicle.N,:).^2;
      % difference of source and target location and distance squared

      rdotfTIMESrdotn = ...
        (diffxy(1:vesicle.N,:).*nx + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).*ny)./dis2.^2 .* ...
        (diffxy(1:vesicle.N,:).*denx + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).*deny);
      % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term
      rdotfTIMESrdotn = rdotfTIMESrdotn * diag(1-vesicle.viscCont(K));
      % have accounted for the scaling with (1-\nu) here

      stokesDLP(j,k) = stokesDLP(j,k) + ...
          sum(sum(rdotfTIMESrdotn.*diffxy(1:vesicle.N,:)));
      stokesDLP(j+vesicle.N,k) = stokesDLP(j+vesicle.N,k) + ...
          sum(sum(rdotfTIMESrdotn.*diffxy(vesicle.N+1:2*vesicle.N,:)));
      % double-layer potential for Stokes
    end
  end

  stokesDLP = stokesDLP/pi;
  % 1/pi is the coefficient in front of the double-layer potential
end
% double-layer potential due to all vesicles except oneself

end % exactStokesDL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLP,stokesDLPtar] = ...
    exactStokesDLfmm(o,vesicle,f,Xtar,K)
% [stokesDLP,stokeDLPtar] = exactStokesDLfmm(vesicle,f,Xtar,K) uses the
% FMM to compute the double-layer potential due to all vesicles except
% itself vesicle is a class of object capsules and f is the density
% function NOT scaled by arclength term.  Xtar is a set of points where
% the double-layer potential due to all vesicles in index set K needs
% to be evaulated
global fmms
disp('OLD FMM is ON')
fmms = fmms + 1;
% count the total number of calls to fmm

oc = curve;
[x,y] = oc.getXY(vesicle.X); % seperate x and y coordinates
nx = vesicle.xt(vesicle.N+1:2*vesicle.N,:);
ny = -vesicle.xt(1:vesicle.N,:);
% seperate the x and y coordinates of the normal vector

den = f.*[vesicle.sa;vesicle.sa]*2*pi/vesicle.N;
den = den * diag(1-vesicle.viscCont);

if (nargin == 5)
  stokesDLP = [];
else
  [f1,f2] = oc.getXY(den);
  % need to multiply by arclength term.  Seperate it into
  % x and y coordinate

  [u,v] = stokesDLPfmm(f1(:),f2(:),x(:),y(:),nx(:),ny(:));

  stokesDLP = zeros(2*vesicle.N,vesicle.nv); % initialize
  for k = 1:vesicle.nv
    is = (k-1)*vesicle.N+1;
    ie = k*vesicle.N;
    stokesDLP(1:vesicle.N,k) = u(is:ie);
    stokesDLP(vesicle.N+1:2*vesicle.N,k) = v(is:ie);
  end
  % Wrap the output of the FMM into the usual 
  % [[x1;y1] [x2;y2] ...] format


  for k = 1:vesicle.nv
    [u,v] = stokesDLPfmm(f1(:,k),f2(:,k),x(:,k),y(:,k),...
        nx(:,k),ny(:,k));
    stokesDLP(:,k) = stokesDLP(:,k) - [u;v];
  end
  % Subtract potential due to each vesicle on its own.  Nothing
  % to change here for potential at Xtar
end

if nargin == 3
  stokesDLPtar = [];
else
  [x,y] = oc.getXY(vesicle.X(:,K)); 
  % seperate x and y coordinates at vesicles indexed by K
  nx = vesicle.xt(vesicle.N+1:2*vesicle.N,K);
  ny = -vesicle.xt(1:vesicle.N,K);
  [Ntar,ncol] = size(Xtar);
  Ntar = Ntar/2;
  x2 = Xtar(1:Ntar,:);
  x = [x(:);x2(:)];
  y2 = Xtar(Ntar+1:2*Ntar,:);
  y = [y(:);y2(:)];
  % Stack the x and y coordinates of the target points
  [f1,f2] = oc.getXY(den(:,K));
  % seperate x and y coordinates at vesicles indexed by K
  f1 = [f1(:);zeros(Ntar*ncol,1)];
  f2 = [f2(:);zeros(Ntar*ncol,1)];
  % pad density function with zeros so that Xtar doesn't
  % affect the double-layer potential
  nx = [nx(:);zeros(Ntar*ncol,1)];
  ny = [ny(:);zeros(Ntar*ncol,1)];
  % pad the normal vector with zeros so that Xtar doesn't
  % affect the double-layer potential

  [u,v] = stokesDLPfmm(f1,f2,x,y,nx,ny);
  
  stokesDLPtar = zeros(2*Ntar,ncol); % initialize
  for k = 1:ncol
    is = vesicle.N*numel(K) + (k-1)*Ntar+1;
    ie = is + Ntar - 1;
    stokesDLPtar(1:Ntar,k) = u(is:ie);
    stokesDLPtar(Ntar+1:2*Ntar,k) = v(is:ie);
  end
  % Wrap the output of the FMM in the usual format
  % for the target points
end

end % exactStokesDLfmm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLP,stokesDLPtar] = ...
    exactStokesDLnewfmm(o,vesicle,f,Xtar,K)

global fmms
fmms = fmms + 1;
% count the total number of calls to fmm

oc = curve;
[x,y] = oc.getXY(vesicle.X); % separate x and y coordinates
nx = vesicle.xt(vesicle.N+1:2*vesicle.N,:);
ny = -vesicle.xt(1:vesicle.N,:);
% separate the x and y coordinates of the normal vector

zn = nx + 1i*ny;
% complex normal vector to the curve

% vector density used in the old fmm
den = f.*[vesicle.sa;vesicle.sa]*2*pi/vesicle.N;
den = den * diag(1-vesicle.viscCont);
denx = den(1:vesicle.N,:);
deny = den(vesicle.N+1:2*vesicle.N,:);

% complex vector density (deny,-denx)
mu = deny - 1i*denx;
dip1 = 1/(4*pi)*mu.*zn;
dip2 = 1/(4*pi)*(mu.*conj(zn)-conj(mu).*zn);

if (nargin == 5)
  stokesDLP = [];
else
  % Get the complex velocity (v,-u)
  vel = stokesDLPnewfmm(dip1(:),dip2(:),x(:),y(:));
  v = real(vel);
  u = -imag(vel);

  stokesDLP = zeros(2*vesicle.N,vesicle.nv); %initialize
  for k = 1:vesicle.nv
    is = (k-1)*vesicle.N+1;
    ie = k*vesicle.N;
    stokesDLP(1:vesicle.N,k) = u(is:ie);
    stokesDLP(vesicle.N+1:2*vesicle.N,k) = v(is:ie);
  end

  for k = 1:vesicle.nv
    vel = stokesDLPnewfmm(dip1(:,k),dip2(:,k),x(:,k),y(:,k));
    v = real(vel); 
    u = -imag(vel);
    stokesDLP(:,k) = stokesDLP(:,k) - [u;v];
  end 
  % subtract the diagonal term  
end

if (nargin == 3)
  stokesDLPtar = [];
else
  [x,y] = oc.getXY(vesicle.X(:,K));
  % separate x and y coordinates of the vesicles indexed by K
  nx = vesicle.xt(vesicle.N+1:2*vesicle.N,K);
  ny = -vesicle.xt(1:vesicle.N,K);
  [Ntar,ncol] = size(Xtar);
  Ntar = Ntar/2;
  x2 = Xtar(1:Ntar,:);
  x = [x(:);x2(:)];
  y2 = Xtar(Ntar+1:2*Ntar,:);
  y = [y(:);y2(:)];
  % Stack the x and y coordinates of the target points
  nx = [nx(:);zeros(Ntar*ncol,1)];
  ny = [ny(:);zeros(Ntar*ncol,1)];
  % pad the normal vector with zeros so that Xtar doesn't
  % affect the double-layer potential

  zn = nx + 1i*ny;
  % complex normal vector

  % vector DLP density
  [denx,deny] = oc.getXY(den(:,K));
  denx = [denx(:);zeros(Ntar*ncol,1)];
  deny = [deny(:);zeros(Ntar*ncol,1)];

  % complex vector DLP density
  mu = deny - 1i*denx;
  dip1 = 1/(4*pi)*mu.*zn;
  dip2 = 1/(4*pi)*(mu.*conj(zn)-conj(mu).*zn);

  vel = stokesDLPnewfmm(dip1(:),dip2(:),x(:),y(:));
  v = real(vel);
  u = -imag(vel);

  stokesDLPtar = zeros(2*Ntar,ncol); % initialize
  for k = 1 : ncol
    is = vesicle.N*numel(K) + (k-1)*Ntar+1;
    ie = is + Ntar-1;
    stokesDLPtar(1:Ntar,k) = u(is:ie);
    stokesDLPtar(Ntar+1:2*Ntar,k) = v(is:ie);
  end
  % Wrap the output of the FMM in the usual format
  % for the target points
end

end % exactStokesDLnewFMM


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF ROUTINES THAT EVALUATE LAYER-POTENTIALS
% WHEN SOURCES ~= TARGETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LP = lagrangeInterp(o)
% interpMap = lagrangeInterp builds the Lagrange interpolation
% matrix that takes seven function values equally distributed
% in [0,1] and returns the seven polynomial coefficients

interpMat = zeros(7);
LP(1,1) = 6.48e1;
LP(1,2) = -3.888e2;
LP(1,3) = 9.72e2;
LP(1,4) = -1.296e3;
LP(1,5) = 9.72e2;
LP(1,6) = -3.888e2;
LP(1,7) = 6.48e1;

LP(2,1) = -2.268e2;
LP(2,2) = 1.296e3;
LP(2,3) = -3.078e3;
LP(2,4) = 3.888e3;
LP(2,5) = -2.754e3;
LP(2,6) = 1.0368e3;
LP(2,7) = -1.62e2;

LP(3,1) = 3.15e2;
LP(3,2) = -1.674e3;
LP(3,3) = 3.699e3;
LP(3,4) = -4.356e3;
LP(3,5) = 2.889e3;
LP(3,6) = -1.026e3;
LP(3,7) = 1.53e2;

LP(4,1) = -2.205e2;
LP(4,2) = 1.044e3;
LP(4,3) = -2.0745e3;
LP(4,4) = 2.232e3;
LP(4,5) = -1.3815e3;
LP(4,6) = 4.68e2;
LP(4,7) = -6.75e1;

LP(5,1) = 8.12e1;
LP(5,2) = -3.132e2;
LP(5,3) = 5.265e2;
LP(5,4) = -5.08e2;
LP(5,5) = 2.97e2;
LP(5,6) = -9.72e1;
LP(5,7) = 1.37e1;

LP(6,1) = -1.47e1;
LP(6,2) = 3.6e1;
LP(6,3) = -4.5e1;
LP(6,4) = 4.0e1;
LP(6,5) = -2.25e1;
LP(6,6) = 7.2e0;
LP(6,7) = -1e0;

LP(7,1) = 1e0;
% rest of the coefficients are zero

end % lagrangeInterp

end % methods 

end % classdef
