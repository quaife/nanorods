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
  profile;
  om;
  
end % properties

methods 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = poten(N, om)
% o = poten(N): constructor; N is the number of points per curve

o.interpMat = o.lagrangeInterp;
% load in the interpolation matrix which is precomputed with 7
% interpolation points
if nargin == 2
    o.om = om;
    o.profile = om.profile;
end

o.N = N;

end % poten: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LP = nearSingInt(o,geomSou,f,selfMat,D,...
    NearStruct,kernel,kernelDirect,geomTar,tEqualS,idebug)
% LP = nearSingInt(geom,f,selfMat,NearStruct,kernel,kernelDirect,
% geomTar,tEqualS,idebug) computes a layer potential due to f at all
% points in geomTar.X.  If tEqualS==true, then the geomTar ==
% geomSou and the self-geom interaction is skipped.  selfMat is
% the diagonal of the potential needed to compute the layer potential of
% each geom indepenedent of all others.  kernel and kernelDirect are
% two (possibly the same) routines that compute the layer potential.
% kernelDirect always uses the direct method whereas kernel may use an
% FMM-accelerated method.  NearStruct is a structure containing the
% variables zone,dist,nearest,icp,argnear which are required by
% near-singular integration (they keep everything sorted and
% precomputed) Everything is in the 2*N x n format Can pass a final
% argument if desired so that plots of the near-singular integration
% algorithm are displayed

if (tEqualS && size(geomSou.X,2) == 1)
  LP = zeros(size(geomSou.X));
  return
end
% only a single geom, so velocity on all other geoms will always
% be zero

if (nargin == 10)
  idebug = false;
end

dist = NearStruct.dist;
zone = NearStruct.zone;
nearest = NearStruct.nearest;
icp = NearStruct.icp;
argnear = NearStruct.argnear;
nearFibers = NearStruct.nearFibers;

Xsou = geomSou.X; % source positions
Nsou = size(Xsou,1)/2; % number of source points
nsou = size(Xsou,2); % number of source 'geoms'
Xtar = geomTar.X; % target positions
Ntar = size(Xtar,1)/2; % number of target points
ntar = size(Xtar,2); % number of target 'geoms'
h = geomSou.length/Nsou; % arclength term

Nup = Nsou*ceil(sqrt(Nsou));

vself = selfMat(f);

% upsample to N^(3/2).  
Xup = [interpft(Xsou(1:Nsou,:),Nup); interpft(Xsou(Nsou+1:2*Nsou,:),Nup)];
fup = [interpft(f(1:Nsou,:),Nup); interpft(f(Nsou+1:2*Nsou,:),Nup)];

geomUp = capsules([],Xup);
% Build an object with the upsampled geom

interpOrder = size(o.interpMat,1);
% lagrange interpolation order
p = ceil((interpOrder+1)/2);
% want half of the lagrange interpolation points to the left of the
% closest point and the other half to the right

if tEqualS % sources == targets
  if nsou > 1
    if (strfind(char(kernel),'fmm'))
      farField = kernel(geomSou,f,D);

      % evaluate layer potential at all targets except ignore the
      % diagonal term
    else
      for k = 1:nsou
        K = [(1:k-1) (k+1:nsou)];
        [~,farField(:,k)] = kernelDirect(geomUp,fup,[],Xtar(:,k),K);
      end
      % This is a huge savings if we are using a direct method rather
      % than the fmm to evaluate the layer potential.  The speedup is
      % more than N^{1/2}, where N is the resolution of the geoms
      % that we are computing with
    end
  else
    farField = zeros(2*Ntar,ntar);
  end

else % sources ~= targets
    [~,farField] = kernel(geomUp,fup,[],Xtar,1:nsou);
    % evaluate layer potential due to all 'geoms' at all points in
    % Xtar
end
% Use upsampled trapezoid rule to compute layer potential

nearField = zeros(2*Ntar,ntar);

beta = 1.1;
% small buffer to make sure Lagrange interpolation points are not in the
% near zone
vel = zeros(2*Ntar, nsou, nsou); %allocate array for vel

for k1 = 1:nsou
    if tEqualS % sources == targets
        % skip diagonal geom
        K = nearFibers{k1};
    else % sources ~= targets
        K = (1:ntar);
        % consider all geoms
    end
    for k2 = K
        J = find(zone{k1}(:,k2) == 1);
        % set of points on geom k2 close to geom k1
        if (numel(J) ~= 0)
            indcp = icp{k1}(J,k2);
            % closest point on geom k1 to each point on geom k2
            % that is close to geom k1
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
            
            if ((numel(J) + numel(f)) >= 512 && numel(J) > 32)
                [~,potTar] = kernel(geomSou,f,D,[Xtar(J,k2);Xtar(J+Ntar,k2)],k1);
            else
                [~,potTar] = kernelDirect(geomSou,f,D,[Xtar(J,k2);Xtar(J+Ntar,k2)],k1);
            end
            % Need to subtract off contribution due to geom k1 since its
            % layer potential will be evaulted using Lagrange interpolant of
            % nearby points
            nearField(J,k2) =  nearField(J,k2) - potTar(1:numel(J));
            nearField(J+Ntar,k2) =  nearField(J+Ntar,k2) - potTar(numel(J)+1:end);
            
            XLag = zeros(2*numel(J),interpOrder - 1);
            % initialize space for initial tracer locations
            for i = 1:numel(J)
                nx = (Xtar(J(i),k2) - nearest{k1}(J(i),k2))/...
                    dist{k1}(J(i),k2);
                ny = (Xtar(J(i)+Ntar,k2) - nearest{k1}(J(i)+Ntar,k2))/...
                    dist{k1}(J(i),k2);
                XLag(i,:) = nearest{k1}(J(i),k2) + ...
                    beta*h(k1)*nx*(1:interpOrder-1);
                XLag(i+numel(J),:) = nearest{k1}(J(i)+Ntar,k2) + ...
                    beta*h(k1)*ny*(1:interpOrder-1);
                % Lagrange interpolation points coming off of geom k1 All
                % points are behind Xtar(J(i),k2) and are sufficiently far from
                % geom k1 so that the Nup-trapezoid rule gives sufficient
                % accuracy
            end
            
            if (numel(XLag)/2 > 100)
                [~,lagrangePts] = kernel(geomUp,fup,[],XLag,k1);
            else
                [~,lagrangePts] = kernelDirect(geomUp,fup,[],XLag,k1);
            end
            % evaluate velocity at the lagrange interpolation points
            
            for i = 1:numel(J)
                Px = o.interpMat*[vel(J(i),k2,k1) ...
                    lagrangePts(i,:)]';
                Py = o.interpMat*[vel(J(i)+Ntar,k2,k1) ...
                    lagrangePts(i+numel(J),:)]';
                % Build polynomial interpolant along the one-dimensional
                % points coming out of the geom
                dscaled = full(dist{k1}(J(i),k2)/(beta*h(k1)*(interpOrder-1)));
                % Point where interpolant needs to be evaluated
                
                v = filter(1,[1 -dscaled],Px);
                nearField(J(i),k2) = nearField(J(i),k2) +  v(end);
                
                v = filter(1,[1 -dscaled],Py);
                nearField(J(i)+Ntar,k2) = nearField(J(i)+Ntar,k2) + v(end);
                
                
                if idebug
                    figure(2); clf; hold on;
                    plot(Xsou(1:Nsou,:),Xsou(Nsou+1:end,:),'r.','markersize',10)
                    plot(Xtar(1:Ntar,:),Xtar(Ntar+1:end,:),'k.','markersize',10)
                    plot(Xtar(J,k2),Xtar(Ntar+J,k2),'b.','markersize',10)
                    plot(XLag(1:numel(J),:),XLag(numel(J)+1:end,:),'kx','markersize',10)
                    plot(XLag(i,:),XLag(numel(J)+i,:),'gx','markersize',10)
                    axis equal
                    
                    figure(1); clf; hold on
                    plot((0:interpOrder-1)*beta*h(k1),...
                        real([vel(J(i),k2,k1) lagrangePts(i,:)]),'g-o')
                    plot((0:interpOrder-1)*beta*h(k1),...
                        real([vel(J(i)+Ntar,k2,k1) lagrangePts(i+numel(J),:)]),'r--o')
                    
                    figure(3)
                    clf
                    hold on
                    plot(f(1:Nsou,k1));
                    plot(f(Nsou+1:2*Nsou,k1));
                    
                    drawnow;
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
% upsampled if source==geom so need to truncate here.  We are 
% only using Ntar target points.  Note that it is only the sources 
% that were upsampled

end % nearSingInt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF ROUTINES THAT BUILD LAYER-POTENTIAL MATRICIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = stokesDLmatrix(~,geom)
% D = stokesDLmatrix(geom), generate double-layer potential for 
% Stokes. geom is a data structure defined as in the capsules class
% D is (2N,2N,n) array where Np is the number of points per curve and 
% n is the number of curves in X 

oc = curve;
[x,y] = oc.getXY(geom.X);
% Vesicle positions
Ng = geom.N;
% number of points per geoms
D = zeros(2*Ng,2*Ng,geom.n);
% initialize space for double-layer potential matrix

for k=1:geom.n  % Loop over curves
  xx = x(:,k);
  yy = y(:,k);
  % locations
  [tx,ty] = oc.getXY(geom.xt);
  tx = tx(:,k); ty = ty(:,k);
  % Vesicle tangent
  sa = geom.sa(:,k)';
  % Jacobian
  cur = geom.cur(:,k)';
  % curvature

  xtar = xx(:,ones(Ng,1));
  ytar = yy(:,ones(Ng,1));
  % target points

  xsou = xx(:,ones(Ng,1))';
  ysou = yy(:,ones(Ng,1))';
  % source points

  txsou = tx';
  tysou = ty';
  % tangent at sources
  sa = sa(ones(Ng,1),:);
  % Jacobian

  diffx = xtar - xsou;
  diffy = ytar - ysou;
  rho4 = (diffx.^2 + diffy.^2).^(-2);
  rho4(1:Ng+1:Ng^2) = 0;
  % set diagonal terms to 0

  kernel = diffx.*(tysou(ones(Ng,1),:)) - ...
            diffy.*(txsou(ones(Ng,1),:));
  kernel = kernel.*rho4.*sa;

  D11 = kernel.*diffx.^2;
  % (1,1) component
  D11(1:Ng+1:Ng^2) = -0.5*cur.*sa(1,:).*txsou.^2;
  % diagonal limiting term

  D12 = kernel.*diffx.*diffy;
  % (1,2) component
  D12(1:Ng+1:Ng^2) = -0.5*cur.*sa(1,:).*txsou.*tysou;
  % diagonal limiting term

  D22 = kernel.*diffy.^2;
  % (2,2) component
  D22(1:Ng+1:Ng^2) = -0.5*cur.*sa(1,:).*tysou.^2;
  % diagonal limiting term

  D(:,:,k) = [D11 D12; D12 D22];
  % build matrix with four blocks
  D(:,:,k) = 1/pi*D(:,:,k)*2*pi/Ng;
  % scale with the arclength spacing and divide by pi
end % k

end % stokesDLmatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N0 = stokesN0matrix(~,wall)
% N0 = stokesN0matrix(vesicle) generates the the integral operator with 
% kernel normal(x) \otimes normal(y) which removes the rank one defficiency 
% of the double-layer potential.  Need this operator for solid walls

normal = [wall.xt(wall.N+1:2*wall.N,:);-wall.xt(1:wall.N,:)]; % Normal vector
normal = normal(:,ones(2*wall.N,1));

sa = [wall.sa(:,1);wall.sa(:,1)];
sa = sa(:,ones(2*wall.N,1));
N0 = zeros(2*wall.N,2*wall.N,wall.n);
N0(:,:,1) = normal.*normal'.*sa'*2*pi/wall.N;
% Use N0 if solving (-1/2 + DLP)\eta = f where f has no flux through
% the boundary.  By solving (-1/2 + DLP + N0)\eta = f, we guarantee
% that \eta also has no flux through the boundary.  This is not
% required, but it means we're enforcing one addition condition on eta
% which removes the rank one kernel.  DLP is the double-layer potential
% for stokes equation

end % stokesN0matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = pressDLmatrix(~,geom)
% P = pressDLmatrix(vesicle), generates the matrix that returns the
% pressure given the traction jump.  Matrix has dimensions (2*N,2*N,n)
% where N is the number of points per curve and n is the number of
% curves in X.  Matrix is not square since traction jump is
% vector-valued whereas the pressure is a scalar-valued function

oc = curve;
[x,y] = oc.getXY(geom.X);
nx = geom.xt(geom.N+1:2*geom.N,:);
ny = -geom.xt(1:geom.N,:);
% Normal vector

P = zeros(geom.N,2*geom.N,geom.n);
% initialize double-layer potential to zero
for k=1:geom.n  % Loop over curves
  index1 = (1:2:geom.N)'; % odd-indexed source points
  for j=2:2:geom.N % Loop over targets
    rho2 = (x(j,k) - x(index1,k)).^2 + (y(j,k) - y(index1,k)).^2;
    % distance squared
    rx = x(j,k) - x(index1,k);
    ry = y(j,k) - y(index1,k);
    rdotn = rx.*nx(index1,k) + ry.*ny(index1,k);

    coeff = (-nx(index1,k)./rho2 + 2*rdotn./rho2.^2.*rx) .* ...
        geom.sa(index1,k)/pi;
    P(j,index1,k) = 4*pi/geom.N*coeff;
    % part that multiplies x-component of traction jump
    % need factor of 4 instead of 2 because we are only using
    % half the terms in the quadrature
    P(j,j,k) = P(j,j,k) - sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2
    coeff = (-ny(index1,k)./rho2 + 2*rdotn./rho2.^2.*ry) .* ...
        geom.sa(index1,k)/pi;
    P(j,geom.N+index1,k) = 4*pi/geom.N*coeff;
    % part that multiplies y-component of traction jump
    % need factor of 4 instead of 2 because we are only using
    % half the terms in the quadrature
    P(j,geom.N+j,k) = P(j,geom.N+j,k) - ...
        sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2
  end % j

  index1 = (2:2:geom.N)'; % even-indexed source points
  for j=1:2:geom.N % Loop over targets
    rho2 = (x(j,k) - x(index1,k)).^2 + (y(j,k) - y(index1,k)).^2;
    % distance squared
    rx = x(j,k) - x(index1,k);
    ry = y(j,k) - y(index1,k);
    rdotn = rx.*nx(index1,k) + ry.*ny(index1,k);
    % dot product of r with normal

    coeff = (-nx(index1,k)./rho2 + 2*rdotn./rho2.^2.*rx) .* ...
        geom.sa(index1,k)/pi;
    P(j,index1,k) = 4*pi/geom.N*coeff;
    % part that multiplies x-component of traction jump
    % need factor of 4 instead of 2 because we are only using
    % half the terms in the quadrature
    P(j,j,k) = P(j,j,k) - sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = (-ny(index1,k)./rho2 + 2*rdotn./rho2.^2.*ry) .* ...
        geom.sa(index1,k)/pi;
    P(j,geom.N+index1,k) = 4*pi/geom.N*coeff;
    % part that multiplies y-component of traction jump
    % need factor of 4 instead of 2 because we are only using
    % half the terms in the quadrature
    P(j,geom.N+j,k) = P(j,geom.N+j,k) - ...
        sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2
  end % j
end % k


end % pressDLmatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D1,D2] = stressDLmatrix(~,geom)
% [D1,D2] = stressDLmatrix(vesicle), generates the matrix that returns
% the stress tensor due to the double-layer potential applied to [1;0] 
% (D1) and to [0;1] (D2) given the traction jump.  Matricies have 
% dimensions (2*N,2*N,n) where N is the number of points per curve 
% and n is the number of curves in X.  Matrix is square since traction 
% jump is vector-valued and so is the stress

normal = [geom.xt(geom.N+1:2*geom.N,:);...
          -geom.xt(1:geom.N,:)]; % Normal vector
oc = curve;
[x,y] = oc.getXY(geom.X);
% Vesicle positions

D1 = zeros(2*geom.N,2*geom.N,geom.n);
D2 = zeros(2*geom.N,2*geom.N,geom.n);
% initialize double-layer potential to zero
for k=1:geom.n  % Loop over curves
  index1 = (1:2:geom.N)'; % odd-indexed source points
  for j=2:2:geom.N % Loop over targets
    rx = x(j,k) - x(index1,k);
    ry = y(j,k) - y(index1,k);
    rho2 = rx.^2 + ry.^2;
    % distance squared
    nx = normal(index1,k);
    ny = normal(index1+geom.N,k);
    % normal vector
    rdotn = rx.*nx + ry.*ny;
    % dot product of r with normal

    coeff = 1./rho2.*nx - ...
        8*rx.*rdotn./rho2.^3.*rx.^2 + ...
        rdotn./rho2.^2.*(2*rx) + ...
        1./rho2.^2.*(2*rx.^2.*nx);
    coeff = coeff/pi.*geom.sa(index1,k);
    D1(j,index1,k) = 4*pi/geom.N*coeff;
    % part of first component of stress1 that multiplies first 
    % component of traction jump
    D1(j,j,k) = D1(j,j,k) - sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = 1./rho2.*ny - ...
        8*rx.*rdotn./rho2.^3.*rx.*ry + ...
        1./rho2.^2.*(2*rx.*ry.*nx);
    coeff = coeff/pi.*geom.sa(index1,k);
    D1(j,index1+geom.N,k) = 4*pi/geom.N*coeff;
    % part of first component of stress1 that multiplies second
    % component of traction jump
    D1(j,geom.N+j,k) = D1(j,geom.N+j,k) - ...
        sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = -8*rx.*rdotn./rho2.^3.*rx.*ry + ...
        rdotn./rho2.^2.*ry + ...
        1./rho2.^2.*(rx.^2.*ny+rx.*ry.*nx);
    coeff = coeff/pi.*geom.sa(index1,k);
    D1(j+geom.N,index1,k) = 4*pi/geom.N*coeff;
    % part of second component of stress1 that multiplies first 
    % component of traction jump
    D1(geom.N+j,j,k) = D1(geom.N+j,j,k) - ...
        sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = -8*rx.*rdotn./rho2.^3.*ry.^2 + ...
        rdotn./rho2.^2.*rx + ...
        1./rho2.^2.*(rx.*ry.*ny+ry.^2.*nx);
    coeff = coeff/pi.*geom.sa(index1,k);
    D1(j+geom.N,index1+geom.N,k) = 4*pi/geom.N*coeff;
    % part of second component of stress1 that multiplies second 
    % component of traction jump
    % Have built the stress tensor applied to [1;0] at half the points
    D1(geom.N+j,geom.N+j,k) = ...
        D1(geom.N+j,geom.N+j,k) - ...
        sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2


    coeff = -8*ry.*rdotn./rho2.^3.*rx.^2 + ...
        rdotn./rho2.^2.*ry + ...
        1./rho2.^2.*(rx.^2.*ny+rx.*ry.*nx);
    coeff = coeff/pi.*geom.sa(index1,k);
    D2(j,index1,k) = 4*pi/geom.N*coeff;
    % part of first component of stress1 that multiplies first 
    % component of traction jump
    D2(j,j,k) = D2(j,j,k) - sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = -8*ry.*rdotn./rho2.^3.*rx.*ry + ...
        rdotn./rho2.^2.*rx + ...
        1./rho2.^2.*(rx.*ry.*ny+ry.^2.*nx);
    coeff = coeff/pi.*geom.sa(index1,k);
    D2(j,index1+geom.N,k) = 4*pi/geom.N*coeff;
    % part of first component of stress1 that multiplies second
    % component of traction jump
    D2(j,geom.N+j,k) = D2(j,geom.N+j,k) - ...
        sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = 1./rho2.*nx - ...
        8*ry.*rdotn./rho2.^3.*rx.*ry + ...
        1./rho2.^2.*(2*rx.*ry.*ny);
    coeff = coeff/pi.*geom.sa(index1,k);
    D2(j+geom.N,index1,k) = 4*pi/geom.N*coeff;
    % part of second component of stress1 that multiplies first 
    % component of traction jump
    D2(geom.N+j,j,k) = D2(geom.N+j,j,k) - ...
        sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = 1./rho2.*ny - ...
        8*ry.*rdotn./rho2.^3.*ry.^2 + ...
        rdotn./rho2.^2.*(2*ry) + ...
        1./rho2.^2.*(2*ry.^2.*ny);
    coeff = coeff/pi.*geom.sa(index1,k);
    D2(j+geom.N,index1+geom.N,k) = 4*pi/geom.N*coeff;
    % part of second component of stress1 that multiplies second 
    % component of traction jump
    D2(geom.N+j,geom.N+j,k) = ...
        D2(geom.N+j,geom.N+j,k) - ...
        sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2
    % Have built the stress tensor applied to [0;1] at half the points


    % need factors of 4 instead of 2 because we are only using
    % half the terms in the quadrature
  end % j

  index1 = (2:2:geom.N)'; % even-indexed source points
  for j=1:2:geom.N % Loop over targets
    rx = x(j,k) - x(index1,k);
    ry = y(j,k) - y(index1,k);
    rho2 = rx.^2 + ry.^2;
    % distance squared
    nx = normal(index1,k);
    ny = normal(index1+geom.N,k);
    % normal vector
    rdotn = rx.*nx + ry.*ny;
    % dot product of r with normal

    coeff = 1./rho2.*nx - ...
        8*rx.*rdotn./rho2.^3.*rx.^2 + ...
        rdotn./rho2.^2.*(2*rx) + ...
        1./rho2.^2.*(2*rx.^2.*nx);
    coeff = coeff/pi.*geom.sa(index1,k);
    D1(j,index1,k) = 4*pi/geom.N*coeff;
    % part of first component of stress1 that multiplies first 
    % component of traction jump
    D1(j,j,k) = D1(j,j,k) - sum(coeff)*4*pi/geom.N;

    coeff = 1./rho2.*ny - ...
        8*rx.*rdotn./rho2.^3.*rx.*ry + ...
        1./rho2.^2.*(2*rx.*ry.*nx);
    coeff = coeff/pi.*geom.sa(index1,k);
    D1(j,index1+geom.N,k) = 4*pi/geom.N*coeff;
    % part of first component of stress1 that multiplies second
    % component of traction jump
    D1(j,geom.N+j,k) = D1(j,geom.N+j,k) - ...
        sum(coeff)*4*pi/geom.N;

    coeff = -8*rx.*rdotn./rho2.^3.*rx.*ry + ...
        rdotn./rho2.^2.*ry + ...
        1./rho2.^2.*(rx.^2.*ny+rx.*ry.*nx);
    coeff = coeff/pi.*geom.sa(index1,k);
    D1(j+geom.N,index1,k) = 4*pi/geom.N*coeff;
    % part of second component of stress1 that multiplies first 
    % component of traction jump
    D1(geom.N+j,j,k) = D1(geom.N+j,j,k) - ...
        sum(coeff)*4*pi/geom.N;

    coeff = -8*rx.*rdotn./rho2.^3.*ry.^2 + ...
        rdotn./rho2.^2.*rx + ...
        1./rho2.^2.*(rx.*ry.*ny+ry.^2.*nx);
    coeff = coeff/pi.*geom.sa(index1,k);
    D1(j+geom.N,index1+geom.N,k) = 4*pi/geom.N*coeff;
    % part of second component of stress1 that multiplies second 
    % component of traction jump
    D1(geom.N+j,geom.N+j,k) = ...
        D1(geom.N+j,geom.N+j,k) - ...
        sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2
    % Have built the stress tensor applied to [0;1] at half the points
    % Have built the stress tensor applied to [1;0] at the other 
    % half of the points


    coeff = -8*ry.*rdotn./rho2.^3.*rx.^2 + ...
        rdotn./rho2.^2.*ry + ...
        1./rho2.^2.*(rx.^2.*ny+rx.*ry.*nx);
    coeff = coeff/pi.*geom.sa(index1,k);
    D2(j,index1,k) = 4*pi/geom.N*coeff;
    % part of first component of stress1 that multiplies first 
    % component of traction jump
    D2(j,j,k) = D2(j,j,k) - sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = -8*ry.*rdotn./rho2.^3.*rx.*ry + ...
        rdotn./rho2.^2.*rx + ...
        1./rho2.^2.*(rx.*ry.*ny+ry.^2.*nx);
    coeff = coeff/pi.*geom.sa(index1,k);
    D2(j,index1+geom.N,k) = 4*pi/geom.N*coeff;
    % part of first component of stress1 that multiplies second
    % component of traction jump
    D2(j,geom.N+j,k) = D2(j,geom.N+j,k) - ...
        sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = 1./rho2.*nx - ...
        8*ry.*rdotn./rho2.^3.*rx.*ry + ...
        1./rho2.^2.*(2*rx.*ry.*ny);
    coeff = coeff/pi.*geom.sa(index1,k);
    D2(j+geom.N,index1,k) = 4*pi/geom.N*coeff;
    % part of second component of stress1 that multiplies first 
    % component of traction jump
    D2(geom.N+j,j,k) = D2(geom.N+j,j,k) - ...
        sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = 1./rho2.*ny - ...
        8*ry.*rdotn./rho2.^3.*ry.^2 + ...
        rdotn./rho2.^2.*(2*ry) + ...
        1./rho2.^2.*(2*ry.^2.*ny);
    coeff = coeff/pi.*geom.sa(index1,k);
    D2(j+geom.N,index1+geom.N,k) = 4*pi/geom.N*coeff;
    % part of second component of stress1 that multiplies second 
    % component of traction jump
    D2(geom.N+j,geom.N+j,k) = ...
        D2(geom.N+j,geom.N+j,k) - ...
        sum(coeff)*4*pi/geom.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2
    % Have built the stress tensor applied to [0;1] at the other 
    % half of the points

    % need factors of 4 instead of 2 because we are only using
    % half the terms in the quadrature
  end % j
end % k

end % stressDLmatrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF ROUTINES THAT BUILD LAYER-POTENTIAL MATRICIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF ROUTINES THAT EVALUATE LAYER-POTENTIALS
% WHEN SOURCES == TARGETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DLP = exactStokesDLdiag(~,geom,D,f)
% DLP = exactStokesDLdiag(geom,f,K) computes the diagonal term of
% the double-layer potential due to f around all geoms.  Source and
% target points are the same.  This uses trapezoid rule with the
% curvature at the diagonal in order to guarantee spectral accuracy.
% This routine can either compute the double-layer potential
% matrix-free, which may upsample the number of source points.  Or, if
% the matrix D is passed in and anti-aliasing is not requested, it will
% simply do the matrix-vector product with the precomputed matrix D.

DLP = zeros(2*geom.N,geom.n);
for k = 1:geom.n
  A = D(:,:,k);
  DLP(:,k) = A * f(:,k);
end
% 
% DLP = permute(sum(bsxfun(@times, D, permute(f,[3 1 2])),2), [1 3 2]);

end % exactStokesDLdiag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N0 = exactStokesN0diag(~,wall,N0,f)
% DLP = exactStokesN0diag(vesicle,f) computes the diagonal term of the
% modification of the double-layer potential due to f around outermost
% wall.  Source and target points are the same.  This uses trapezoid
% rule
if isempty(N0)
  Nf = size(f,1)/2;
  oc = curve;
  [fx,fy] = oc.getXY(f(:,1));
  fx = fx.*wall.sa(:,1);
  fy = fy.*wall.sa(:,1);
  [tx,ty] = oc.getXY(wall.xt(:,1));
  % tangent vector
  const = sum(ty.*fx - tx.*fy)*2*pi/Nf;
  % function to be integrated is dot product of normal with density
  % function
  N0 = zeros(2*Nf,1);
  N0 = const*[ty;-tx];
else
  N0 = N0(:,:,1)*f(:,1);
end

end % exactStokesN0diag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pressure = exactPressureDLdiag(~,geom,P,f)
% pressure = exactPressureDLdiag(vesicle,P,f) computes the diagonal
% term of the pressure of the double-layer potental due to f around
% each vesicle.  Source and target points are the same.  For now, we
% just pass the matrix for the layer potential and loop over the
% vesicles

pressure = zeros(size(f));
for k = 1:geom.n
  pressure(:,k) = P(:,:,k) * f(:,k);
end

end % exactPressureDLdiag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stress = exactStressDLdiag(~,geom,S,f)
% stress = exactStressSLdiag(vesicle,S,f) computes the diagonal term of
% the stress of the double-layer potental due to f around each
% vesicle.  Source and target points are the same.  For now, we just
% pass the matrix for the layer potential and loop over the vesicles

stress = zeros(size(f));
for k = 1:geom.n
  stress(:,k) = S(:,:,k) * f(:,k);
end

end % exactStressDLdiag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF ROUTINES THAT EVALUATE LAYER-POTENTIALS
% WHEN SOURCES == TARGETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF ROUTINES THAT EVALUATE LAYER-POTENTIALS
% WHEN SOURCES ~= TARGETS.  CAN COMPUTE LAYER POTENTIAL ON EACH
% VESICLE DUE TO ALL OTHER VESICLES (ex. stokesSLP) AND CAN
% COMPUTE LAYER POTENTIAL DUE TO VESICLES INDEXED IN K1 AT 
% TARGET POINTS Xtar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLP,stokesDLPtar] = exactStokesDL(~, geom, f, ~, Xtar,K1)
% [stokesDLP,stokesDLPtar] = exactStokesDL(geom,f,Xtar,K1) computes the
% double-layer potential due to f around all parts of the geometry
% except itself.  Also can pass a set of target points Xtar and a
% collection of geom K1 and the double-layer potential due to components
% of the geometry in K1 will be evaluated at Xtar.  Everything but Xtar
% is in the 2*N x n format Xtar is in the 2*Ntar x ncol format

if nargin == 6
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
den = f.*[geom.sa;geom.sa]*2*pi/geom.N;
% jacobian term and 2*pi/N accounted for here

oc = curve;
[xsou,ysou] = oc.getXY(geom.X(:,K1));
xsou = xsou(:); ysou = ysou(:);
xsou = xsou(:,ones(Ntar,1))';
ysou = ysou(:,ones(Ntar,1))';

[denx,deny] = oc.getXY(den(:,K1));
denx = denx(:); deny = deny(:);
denx = denx(:,ones(Ntar,1))';
deny = deny(:,ones(Ntar,1))';

[normalx,normaly] = oc.getXYperp(geom.xt(:,K1));
normalx = normalx(:); normaly = normaly(:);
normalx = normalx(:,ones(Ntar,1))';
normaly = normaly(:,ones(Ntar,1))';

for k = 1:ncol % loop over columns of target points
  [xtar,ytar] = oc.getXY(Xtar(:,k));
  xtar = xtar(:,ones(geom.N*numel(K1),1));
  ytar = ytar(:,ones(geom.N*numel(K1),1));
  
  diffx = xtar-xsou; diffy = ytar-ysou;
  dis2 = (diffx).^2 + (diffy).^2;
  % difference of source and target location and distance squared
  
  rdotnTIMESrdotf = (diffx.*normalx + diffy.*normaly)./dis2.^2 .* ...
      (diffx.*denx + diffy.*deny);
  % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term
  
  stokesDLPtar(1:Ntar,k) = stokesDLPtar(1:Ntar,k) + ...
      sum(rdotnTIMESrdotf.*diffx,2);
  stokesDLPtar(Ntar+1:end,k) = stokesDLPtar(Ntar+1:end,k) + ...
      sum(rdotnTIMESrdotf.*diffy,2);
  % r \otimes r term of the double-layer potential
end
stokesDLPtar = stokesDLPtar/pi;
% double-layer potential due to geometry components indexed over K1
% evaluated at arbitrary points

stokesDLP = zeros(2*geom.N,geom.n);
if (nargin == 4 && geom.n > 1)
    for k = 1:geom.n
        K = [(1:k-1) (k+1:geom.n)];
        [x,y] = oc.getXY(geom.X(:,K));
        [nx,ny] = oc.getXYperp(geom.xt(:,K));
        [denx,deny] = oc.getXY(den(:,K));
        for j=1:geom.N
            diffxy = [geom.X(j,k) - x ; geom.X(j+geom.N,k) - y];
            dis2 = diffxy(1:geom.N,:).^2 + ...
                diffxy(geom.N+1:2*geom.N,:).^2;
            % difference of source and target location and distance squared
            
            rdotfTIMESrdotn = ...
                (diffxy(1:geom.N,:).*nx + ...
                diffxy(geom.N+1:2*geom.N,:).*ny)./dis2.^2 .* ...
                (diffxy(1:geom.N,:).*denx + ...
                diffxy(geom.N+1:2*geom.N,:).*deny);
            % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term
            
            stokesDLP(j,k) = stokesDLP(j,k) + ...
                sum(sum(rdotfTIMESrdotn.*diffxy(1:geom.N,:)));
            stokesDLP(j+geom.N,k) = stokesDLP(j+geom.N,k) + ...
                sum(sum(rdotfTIMESrdotn.*diffxy(geom.N+1:2*geom.N,:)));
            % double-layer potential for Stokes
        end
    end
    
    stokesDLP = stokesDLP/pi;
    % 1/pi is the coefficient in front of the double-layer potential
end
% double-layer potential due to all components of the geometry except
% oneself

end % exactStokesDL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLP,stokesDLPtar] = exactStokesDLfmm(o,geom,f, D, Xtar,K)
    % [stokesDLP,stokeDLPtar] = exactStokesDLfmm(geom,f,Xtar,K) uses the
    % FMM to compute the double-layer potential due to all geoms except
    % itself geom is a class of object capsules and f is the density
    % function NOT scaled by arclength term.  Xtar is a set of points where
    % the double-layer potential due to all geoms in index set K needs
    % to be evaulated
    
    oc = curve;
%     [x,y] = oc.getXY(geom.X); % seperate x and y coordinates
%     nx = geom.xt(geom.N+1:2*geom.N,:);
%     ny = -geom.xt(1:geom.N,:);
    % seperate the x and y coordinates of the normal vector
    
    %f = ones(size(f));
    den = f.*[geom.sa;geom.sa]*2*pi/geom.N;
    
    if (nargin == 6)
        stokesDLP = [];
    else
        % upsample source to N^(3/2).
        Nup = geom.N*ceil(sqrt(geom.N));
        Xtar = geom.X;
        
        Xup= [interpft(Xtar(1:geom.N,:),Nup); interpft(Xtar(geom.N+1:end,:),Nup)];
        fup = [interpft(f(1:geom.N,:),Nup); interpft(f(geom.N+1:end,:),Nup)];
        geomUp = capsules([],Xup);
        
        denUp = fup.*[geomUp.sa;geomUp.sa]*2*pi/geomUp.N;
        
        [x,y] = oc.getXY(Xup); % seperate x and y coordinates
    
        nx = geomUp.xt(geomUp.N+1:end,:);
        ny = -geomUp.xt(1:geomUp.N,:);
    
        stokesDLPup = zeros(2*geomUp.N,geomUp.n);
        [fx,fy] = oc.getXY(denUp);
        
        dip1 = 0.25/pi*(fy - 1i*fx).*(nx + 1i*ny);
        dip2 = -1i*0.5/pi*(fx.*nx + fy.*ny);
        
        if o.profile
            tfmmTotal = tic;
        end
        
        % FMM on upsampled data
        velUp = stokesDLPfmm(dip1(:),dip2(:),x(:),y(:));
        u = -imag(velUp);
        v = real(velUp);
        
        if o.profile
            o.om.writeMessage(['FMM for all points took ', ...
                num2str(toc(tfmmTotal)), ' seconds']);
        end
        
        % Wrap the output of the FMM into the usual
        % [[x1;y1] [x2;y2] ...] format
        for k = 1:geom.n
            is = (k-1)*geomUp.N+1;
            ie = k*geomUp.N;
            stokesDLPup(1:geomUp.N,k) = u(is:ie);
            stokesDLPup(geomUp.N+1:end,k) = v(is:ie);
        end
        
        [tx,ty] = oc.getXY(geomUp.xt);

        % ADD DIAGONAL TERM
        for k = 1:geom.n %in principle this could be done using a parfor loop
            txk = tx(:,k); tyk = ty(:,k);
            sa = geomUp.sa(:,k);
            cur = geomUp.cur(:,k);
            
            fDotTau = txk.*fup(1:end/2,k) + tyk.*fup(end/2+1:end,k);
            diag = -[fDotTau.*cur.*sa.*txk; fDotTau.*cur.*sa.*tyk]/geomUp.N;
            
            stokesDLPup(:,k) = stokesDLPup(:,k) + diag;
        end
        
        % downsample to original N points
        stokesDLP = stokesDLPup(1:geomUp.N/geom.N:end,:);
%         stokesDLP =  [interpft(stokesDLPup(1:geomUp.N,:),geom.N); ...
%                     interpft(stokesDLPup(geomUp.N+1:end,:),geom.N)]; 
                
        % Subtract potential due to each geom on its own.  Nothing
        % to change here for potential at Xtar        
        diagDL = o.exactStokesDLdiag(geom, D, f);        
        
        if o.profile
            tfmmLoop = tic;
        end
        
        for k = 1:geom.n %in principle this could be done using a parfor loop
            
            stokesDLP(:,k) = stokesDLP(:,k) - diagDL(:,k);
        end
        
        if o.profile
            o.om.writeMessage(['Direct summation for self interactions took ', ...
                num2str(toc(tfmmLoop)), ' seconds']);
        end
    end
    
    if nargin == 4
        stokesDLPtar = [];
    else

        [Ntar,ncol] = size(Xtar);
        Ntar = Ntar/2;
        stokesDLPtar = zeros(2*Ntar,ncol); % initialize
        [x,y] = oc.getXY(geom.X(:,K));
        % seperate x and y coordinates at geoms indexed by K
        nx = geom.xt(geom.N+1:2*geom.N,K);
        ny = -geom.xt(1:geom.N,K);
        x2 = Xtar(1:Ntar,:);
        x = [x(:);x2(:)];
        y2 = Xtar(Ntar+1:2*Ntar,:);
        y = [y(:);y2(:)];
        % Stack the x and y coordinates of the target points
        [fx,fy] = oc.getXY(den(:,K));
        % seperate x and y coordinates at geoms indexed by K
        
        fx = [fx(:);zeros(Ntar*ncol,1)];
        fy = [fy(:);zeros(Ntar*ncol,1)];
        % pad density function with zeros so that Xtar doesn't
        % affect the double-layer potential
        nx = [nx(:);zeros(Ntar*ncol,1)];
        ny = [ny(:);zeros(Ntar*ncol,1)];
        % pad the normal vector with zeros so that Xtar doesn't
        % affect the double-layer potential
        
        dip1 = 0.25/pi*(fy - 1i*fx).*(nx + 1i*ny);
        dip2 = -1i*0.5/pi*(fx.*nx + fy.*ny);
        
        vel = stokesDLPfmm(dip1(:),dip2(:),x(:),y(:));
        u = -imag(vel);
        v = real(vel);
        
        for k = 1:ncol
            is = geom.N*numel(K) + (k-1)*Ntar+1;
            ie = is + Ntar - 1;
            stokesDLPtar(1:Ntar,k) = u(is:ie);
            stokesDLPtar(Ntar+1:2*Ntar,k) = v(is:ie);
        end
        % Wrap the output of the FMM in the usual format
        % for the target points
    end
    
end % exactStokesDLfmm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [laplaceDLP,laplaceDLPtar] =...
                        exactLaplaceDL(~,geom,f, ~, Xtar,K1)
% pot = exactLaplaceDL(vesicle,f,Xtar,K1) computes the double-layer
% laplace potential due to f around all vesicles except itself.  Also
% can pass a set of target points Xtar and a collection of vesicles K1
% and the double-layer potential due to vesicles in K1 will be
% evaluated at Xtar.  Everything but Xtar is in the 2*N x n format
% Xtar is in the 2*Ntar x ncol format

oc = curve;

nx = geom.xt(geom.N+1:2*geom.N,:);
ny = -geom.xt(1:geom.N,:);

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  laplaceDLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  laplaceDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, user does not need the layer potential at arbitrary
  % points
end

den = f.*[geom.sa;geom.sa]*2*pi/geom.N;
% multiply by arclength term

[xsou,ysou] = oc.getXY(geom.X(:,K1));
xsou = xsou(:); ysou = ysou(:);
xsou = xsou(:,ones(Ntar,1))';
ysou = ysou(:,ones(Ntar,1))';

[denx,deny] = oc.getXY(den(:,K1));
denx = denx(:); deny = deny(:);
denx = denx(:,ones(Ntar,1))';
deny = deny(:,ones(Ntar,1))';

nxK1 = nx(:,K1); nyK1 = ny(:,K1);
nxK1 = nxK1(:); nyK1 = nyK1(:);
nxK1 = nxK1(:,ones(Ntar,1))';
nyK1 = nyK1(:,ones(Ntar,1))';

for k2 = 1:ncol % loop over columns of target points
  [xtar,ytar] = oc.getXY(Xtar(:,k2));
  xtar = xtar(:,ones(geom.N*numel(K1),1));
  ytar = ytar(:,ones(geom.N*numel(K1),1));
  
  diffx = xsou-xtar; diffy = ysou-ytar;
  dis2 = diffx.^2 + diffy.^2;
  
  coeff = (diffx.*nxK1 + diffy.*nyK1)./dis2;
  
  val = coeff.*denx;
  laplaceDLPtar(1:Ntar,k2) = sum(val,2);
  
  val = coeff.*deny;
  laplaceDLPtar(Ntar+1:2*Ntar,k2) = sum(val,2);
end % end k2
% Evaluate double-layer potential at arbitrary target points
laplaceDLPtar = 1/(2*pi)*laplaceDLPtar;
% 1/2/pi is the coefficient in front of the double-layer potential

laplaceDLP = zeros(2*geom.N,geom.n); % Initialize to zero
% if we only have one vesicle, vesicles of course can not collide
% Don't need to run this loop in this case
if (nargin == 3 && geom.n > 1)
  for k1 = 1:geom.n % number of bodies
    K = [(1:k1-1) (k1+1:geom.n)];
    % Loop over all vesicles except k1

    [xsou,ysou] = oc.getXY(geom.X(:,K));
    xsou = xsou(:); ysou = ysou(:);
    xsou = xsou(:,ones(geom.N,1))';
    ysou = ysou(:,ones(geom.N,1))';
    
    [denxK,denyK] = oc.getXY(den(:,K));
    denxK = denxK(:); 
    denxK = denxK(:,ones(geom.N,1))';
    denyK = denyK(:); 
    denyK = denyK(:,ones(geom.N,1))';
    
    nxK = nx(:,K); nyK = ny(:,K);
    nxK = nxK(:); nyK = nyK(:);
    nxK = nxK(:,ones(geom.N,1))';
    nyK = nyK(:,ones(geom.N,1))';
    
    [xtar,ytar] = oc.getXY(geom.X(:,k1));
    xtar = xtar(:); ytar = ytar(:);
    xtar = xtar(:,ones(geom.N*numel(K),1));
    ytar = ytar(:,ones(geom.N*numel(K),1));
    
    diffx = xsou-xtar; diffy = ysou-ytar;
    dis2 = diffx.^2 + diffy.^2;
    
    coeff = (diffx.*nxK + diffy.*nyK)./dis2;

    val = coeff.*denxK;
    laplaceDLP(1:geom.N,k1) = sum(val,2);
    % Laplace DLP with x-coordinate of density function
    val = coeff.*denyK;
    laplaceDLP(geom.N+1:end,k1) = sum(val,2);
    % Laplace DLP with y-coordinate of density function
  end % k1
  % Evaluate double-layer potential at vesicles but oneself
  laplaceDLP = 1/(2*pi)*laplaceDLP;
  % 1/2/pi is the coefficient in front of the double-layer potential
end % nargin == 3

end % exactLaplaceDL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [laplaceDLP,laplaceDLPtar] = exactLaplaceDLfmm(~,geom,f,~,Xtar,K)
% [laplaceDLP,laplaceDLPtar] = exactLaplaceDLfmm(vesicle,f,Xtar,K) uses
% the FMM to compute the double-layer potential due to all vesicles
% except itself vesicle is a class of object capsules and f is the
% density function NOT scaled by arclength term.  Xtar is a set of
% points where the double-layer potential due to all vesicles in index
% set K needs to be evaulated

oc = curve;

den = f.*[geom.sa;geom.sa]*2*pi/geom.N;

if (nargin == 6)
  laplaceDLP = [];
else
  [x,y] = oc.getXY(geom.X); % seperate x and y coordinates
  [tx,ty] = oc.getXY(geom.xt); % tangent vector
  nx = ty; ny = -tx; % normal vector
  [f1,f2] = oc.getXY(den);
  % need to multiply by arclength term.  Seperate it into
  % x and y coordinate

  potx = laplaceDLPfmm(f1(:),x(:),y(:),nx(:),ny(:));
  poty = laplaceDLPfmm(f2(:),x(:),y(:),nx(:),ny(:));
  laplaceDLP = zeros(2*geom.N,geom.n); % initialize
  for k = 1:geom.n
    is = (k-1)*geom.N+1;
    ie = k*geom.N;
    laplaceDLP(1:geom.N,k) = potx(is:ie);
    laplaceDLP(geom.N+1:2*geom.N,k) = poty(is:ie);
  end
  % Wrap the output of the FMM into the usual 
  % [[x1;y1] [x2;y2] ...] format

  for k = 1:geom.n
    potx = laplaceDLPfmm(f1(:,k),x(:,k),y(:,k),nx(:,k),ny(:,k));
    poty = laplaceDLPfmm(f2(:,k),x(:,k),y(:,k),nx(:,k),ny(:,k));
    laplaceDLP(:,k) = laplaceDLP(:,k) - [potx;poty];
  end
  % Subtract potential due to each vesicle on its own.  Nothing
  % to change here for potential at Xtar
end

if nargin == 4
  laplaceDLPtar = [];
else
  [x,y] = oc.getXY(geom.X(:,K)); % seperate x and y coordinates
  [tx,ty] = oc.getXY(geom.xt(:,K)); % tangent vector
  nx = ty; ny = -tx; % normal vector
  % seperate x and y coordinates at vesicles indexed by K
  [Ntar,ncol] = size(Xtar);
  Ntar = Ntar/2;
  x2 = Xtar(1:Ntar,:);
  x = [x(:);x2(:)];
  y2 = Xtar(Ntar+1:2*Ntar,:);
  y = [y(:);y2(:)];
  % Stack the x and y coordinates of the target points
  [f1,f2] = oc.getXY(den(:,K));
  % seperate x and y coordinates at vesicles indexed by K
  nx = [nx(:);zeros(Ntar*ncol,1)];
  ny = [ny(:);zeros(Ntar*ncol,1)];
  f1 = [f1(:);zeros(Ntar*ncol,1)];
  f2 = [f2(:);zeros(Ntar*ncol,1)];
  % pad density function with zeros so that Xtar doesn't
  % affect the single-layer potential

  potx = laplaceDLPfmm(f1,x,y,nx,ny);
  poty = laplaceDLPfmm(f2,x,y,nx,ny);
  laplaceDLPtar = zeros(2*Ntar,ncol); % initialize
  for k = 1:ncol
    is = geom.N*numel(K) + (k-1)*Ntar+1;
    ie = is + Ntar - 1;
    laplaceDLPtar(1:Ntar,k) = potx(is:ie);
    laplaceDLPtar(Ntar+1:2*Ntar,k) = poty(is:ie);
  end
  % Wrap the output of the FMM in the usual format
  % for the target points
end

end % exactLaplaceDLfmm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF ROUTINES THAT EVALUATE LAYER-POTENTIALS
% WHEN SOURCES ~= TARGETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pressDLP,pressDLPtar] = exactPressDL(~,geom,f,~,Xtar,K1)
% [pressDLP,pressDLPtar] = exactPressDL(vesicle,f,pressTrap,Xtar,K1)
% computes the pressure due to all vesicles contained in vesicle and
% indexed over K1.  Evaluates it at Xtar Everything but Xtar is in the
% 2*N x n format Xtar is in the 2*Ntar x ncol format

den = f.*[geom.sa;geom.sa]*2*pi/geom.N;
oc = curve;
[x,y] = oc.getXY(geom.X);
[denx,deny] = oc.getXY(den);
nx = geom.xt(geom.N+1:2*geom.N,:);
ny = -geom.xt(1:geom.N,:);


if nargin == 6
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  pressDLPtar = zeros(Ntar,ncol);
else
  K1 = [];
  pressDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 6, user does not need the layer potential at arbitrary
  % points
end

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis2 = (Xtar(j,k2) - x(:,K1)).^2 + (Xtar(j+Ntar,k2) - y(:,K1)).^2;
    diffxy = [Xtar(j,k2) - x(:,K1) ; Xtar(j+Ntar,k2) - y(:,K1)];
    % distance squared and difference of source and target location
    rdotn = diffxy(1:geom.N,:).*nx(:,K1) + ...
        diffxy(geom.N+1:2*geom.N,:).*ny(1:geom.N,K1); 
  
    val = (nx(:,K1) - 2*rdotn./dis2.*...
                diffxy(1:geom.N,:))./dis2 .* denx(:,K1);
    val = val + (ny(:,K1) - 2*rdotn./dis2.*...
                 diffxy(geom.N+1:2*geom.N,:))./dis2 .* deny(:,K1);

    pressDLPtar(j,k2) = sum(val(:)); 
  end % j
end % k2
% pressure coming from the double-layer potential for Stokes flow

pressDLP = zeros(geom.N,geom.n);
% TODO: NOT SURE IF WE WILL EVER NEED THIS BUT SHOULD PUT IT
% IN NONETHELESS

pressDLPtar = -[pressDLPtar;pressDLPtar]*1/pi;
% near-singular integration needs vector-valued functions also need to
% multiply by 1/(2*pi) as per the pressure of the single-layer
% potential

end % exactPressDL


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stressDLP,stressDLPtar] = exactStressDL1(~,geom,f,~,Xtar,K1)
% [stressDLP,stressDLPtar] = exactStressDL1(vesicle,f,Xtar,K1) computes
% the stress due to the double-layer potential of all vesicles
% contained in vesicle and indexed over K1.  Only computes the stress
% applied to the direction [1;0].  Evaluates it at Xtar. Everything but
% Xtar is in the 2*N x n format Xtar is in the 2*Ntar x ncol format

oc = curve;
[x,y] = oc.getXY(geom.X);
% Vesicle positions
nx = geom.xt(geom.N+1:2*geom.N,:);
ny = -geom.xt(1:geom.N,:);
% normal componenets

if nargin == 6
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stressDLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stressDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 6, user does not need the layer potential at arbitrary points
end

[fx,fy] = oc.getXY(f);
% first and second components of the density function

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis2 = (Xtar(j,k2) - x(:,K1)).^2 + (Xtar(j+Ntar,k2) - y(:,K1)).^2;
    diffxy = [Xtar(j,k2) - x(:,K1) ; Xtar(j+Ntar,k2) - y(:,K1)];
    % distance squared and difference of source and target location
    [rx,ry] = oc.getXY(diffxy);

    rdotf = rx.*fx(:,K1) + ry.*fy(:,K1);
    % dot product of r and f
    fdotn = fx(:,K1).*nx(:,K1) + fy(:,K1).*ny(:,K1);
    % dot product of f and n
    rdotn = rx.*nx(:,K1) + ry.*ny(:,K1);
    % dot product of r and n

    val = (fdotn./dis2 - ...
        8./dis2.^3.*rdotn.*rdotf.*rx.*rx + ...
        rdotn./dis2.^2.*(2*rx.*fx(:,K1)) + ...
        rdotf./dis2.^2.*(2*rx.*nx(:,K1))).*geom.sa(:,K1);
    % first component of the stress of the double-layer potential
    % applied to [1;0]

    stressDLPtar(j,k2) = sum(val(:)); 
    % scale by arclength
    stressDLPtar(j,k2) = stressDLPtar(j,k2)*2*pi/geom.N;
    % d\theta term

    val = (-8./dis2.^3.*rdotn.*rdotf.*rx.*ry + ...
        rdotn./dis2.^2.*(rx.*fy(:,K1) + ry.*fx(:,K1)) + ...
        rdotf./dis2.^2.*(rx.*ny(:,K1) + ry.*nx(:,K1))).*...
        geom.sa(:,K1);
    % second component of the stress of the double-layer potential
    % applied to [1;0]

    stressDLPtar(j+Ntar,k2) = sum(val(:)); 
    % scale by arclength
    stressDLPtar(j+Ntar,k2) = stressDLPtar(j+Ntar,k2)*2*pi/geom.N;
    % d\theta term

  end % j
end % k2
% stress coming from the single-layer potential for Stokes flow

stressDLP = zeros(2*geom.N,geom.n);
% TODO: NOT SURE IF WE WILL EVER NEED THIS BUT SHOULD PUT IT
% IN NONETHELESS

stressDLPtar = stressDLPtar/pi;
% 1/pi is the constant in front of the stress of the double-layer 
% potential

end % exactStressDL1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stressDLP,stressDLPtar] = exactStressDL2(~,geom,f,~,Xtar,K1)
% [stressDLP,stressDLPtar] = exactStressDL2(vesicle,f,Xtar,K1) computes
% the stress due to the double-layer potential of all vesicles
% contained in vesicle and indexed over K1.  Only computes the stress
% applied to the direction [0;1].  Evaluates it at Xtar Everything but
% Xtar is in the 2*N x n format Xtar is in the 2*Ntar x ncol format

oc = curve;
[x,y] = oc.getXY(geom.X);
% Vesicle positions
nx = geom.xt(geom.N+1:2*geom.N,:);
ny = -geom.xt(1:geom.N,:);
% normal componenets

if nargin == 6
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stressDLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stressDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 6, user does not need the layer potential at arbitrary points
end

[fx,fy] = oc.getXY(f);
% first and second components of the density function

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis2 = (Xtar(j,k2) - x(:,K1)).^2 + (Xtar(j+Ntar,k2) - y(:,K1)).^2;
    diffxy = [Xtar(j,k2) - x(:,K1) ; Xtar(j+Ntar,k2) - y(:,K1)];
    % distance squared and difference of source and target location
    [rx,ry] = oc.getXY(diffxy);

    rdotf = rx.*fx(:,K1) + ry.*fy(:,K1);
    % dot product of r and f
    fdotn = fx(:,K1).*nx(:,K1) + fy(:,K1).*ny(:,K1);
    % dot product of f and n
    rdotn = rx.*nx(:,K1) + ry.*ny(:,K1);
    % dot product of r and n

    val = (-8./dis2.^3.*rdotn.*rdotf.*ry.*rx + ...
        rdotn./dis2.^2.*(rx.*fy(:,K1) + ry.*fx(:,K1)) + ...
        rdotf./dis2.^2.*(rx.*ny(:,K1) + ry.*nx(:,K1))).*...
        geom.sa(:,K1);
    % second component of the stress of the double-layer potential
    % applied to [0;1]

    stressDLPtar(j,k2) = sum(val(:)); 
    stressDLPtar(j,k2) = stressDLPtar(j,k2)*2*pi/geom.N;
    % d\theta term

    val = (fdotn./dis2 - ...
        8./dis2.^3.*rdotn.*rdotf.*ry.*ry + ...
        rdotn./dis2.^2.*(2*ry.*fy(:,K1)) + ...
        rdotf./dis2.^2.*(2*ry.*ny(:,K1))).*...
        geom.sa(:,K1);
    % first component of the stress of the double-layer potential
    % applied to [0;1]

    stressDLPtar(j+Ntar,k2) = sum(val(:)); 
    stressDLPtar(j+Ntar,k2) = stressDLPtar(j+Ntar,k2)*2*pi/geom.N;
    % d\theta term
  end % j
end % k2
% stress coming from the single-layer potential for Stokes flow

stressDLP = zeros(2*geom.N,geom.n);
% TODO: NOT SURE IF WE WILL EVER NEED THIS BUT SHOULD PUT IT
% IN NONETHELESS

stressDLPtar = stressDLPtar/pi;
% 1/pi is the constant in front of the stress of the double-layer 
% potential

end % exactStressDL2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % methods
methods(Static)

    
function LP = lagrangeInterp(~)
% interpMap = lagrangeInterp builds the Lagrange interpolation
% matrix that takes seven function values equally distributed
% in [0,1] and returns the seven polynomial coefficients

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

end % methods(Static)

end % classdef
