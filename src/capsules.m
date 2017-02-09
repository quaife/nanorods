classdef capsules < handle
% This class implements standard calculations that need to be done to a
% componenet of the solid wall, or a collection of arbitrary target
% points (such as tracers).  The main tasks that can be performed
% include constructing structures required for near-singluar integration

properties
N;        % number of points per componenet
nv;       % number of components
X;        % positions of component
center;   % center of each rigid body
xt;       % tangent unit vector
sa;       % Jacobian
cur;      % curvature
length;   % total length of each component
end %properties


methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = capsules(prams, varargin)
% capsules(X) constructs an object of class capsules.  Mostly, it
% takes the values of prams and options that it requires.
% This is the constructor

oc = curve;

if length(varargin) == 1 % coordinates provided directly
    o.X = varargin{1};
    o.N = size(o.X,1)/2; % points per component
    o.nv = size(o.X,2);  % number of components
    o.center = [mean(o.X(1:end/2,:));mean(o.X(end/2+1:end,:))];
    % centers of each component
    
else if length(varargin) == 2 % centres and orientation angles
        
        if prams.nv > 0
            xc = varargin{1};
            tau = varargin{2};
            
            o.N = prams.N;              % points per component
            o.nv = prams.nv;            % number of components
            o.center = xc;
            
            theta = (0:prams.N-1)'*2*pi/prams.N;
            o.X = zeros(2*prams.N,prams.nv);
            
            r = (cos(-theta).^prams.order + sin(-theta).^prams.order).^(-1/prams.order);
            
            for k = 1:prams.nv
                
                x_square = prams.lengths(1)/2*r.*cos(-theta);
                y_square = prams.widths(1)/2*r.*sin(-theta);
                
                o.X(:,k) = [x_square*cos(tau(k)) - y_square*sin(tau(k)) + xc(1,k);
                    x_square*sin(tau(k)) + y_square*cos(tau(k)) + xc(2,k)];
                
                for i = 1:2
                    o.X(:,k) = oc.redistributeArcLength(o.X(:,k));
                end
            end
        else
            o.X = [];
        end
        
    end
end

if ~isempty(o.X)
    [o.sa,o.xt,o.cur] = oc.diffProp(o.X);
    [~,o.length] = oc.geomProp(o.X);
else
    o.sa = [];
    o.xt = [];
    o.cur = [];
    o.length = [];
    o.N = 0;
    o.nv = 0;
    o.center = [];
end

end % capsules: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew, xc, tau, add] = add_fibre(o, rmin, rmax, prams, k)
    max_attempts = 1000;
        
    attempt = 0;
    
    theta = (0:prams.N-1)'*2*pi/prams.N;    
    Xoriginal = o.X;
    rt = (cos(theta).^prams.order+ sin(theta).^prams.order).^(-1/prams.order);
    
    x_square = mean(prams.lengths)/2*rt.*cos(theta);
    y_square = mean(prams.widths)/2*rt.*sin(theta);

    while attempt < max_attempts
        
        add = true;
        
        r = rmin + (rmax-rmin)*rand(1,1);
        theta = 2*pi*rand(1,1);
        
        xc = [r*cos(theta);r*sin(theta)];
        tau = 2*pi*rand(1,1);
        
        o.X = [Xoriginal, [x_square*cos(tau) - y_square*sin(tau) + xc(1);
            x_square*sin(tau) + y_square*cos(tau) + xc(2)]];
        o.nv = o.nv + 1;
        
        [near,~] = o.getZone(o,1);
        
        rx = sqrt(o.X(1:end/2,end).^2+o.X(end/2+1:end,end).^2);
        
        for i = 1:o.nv
            
            if (length(near.nearFibers{i}) > 0 || min(rx) < rmin + 0.01 || max(rx) > rmax - 0.01)
                add = false;                
            end            
        end
        
        attempt = attempt + 1;
        o.X = Xoriginal;
        o.nv = o.nv - 1;
        
        if add
            Xnew = [x_square*cos(tau) - y_square*sin(tau) + xc(1);
                    x_square*sin(tau) + y_square*cos(tau) + xc(2)];
            
            disp(['inserted fibre ', num2str(k)]);
            break;
        else
            disp(['failed to insert fibre ', num2str(k), '(attempt ', num2str(attempt), '), trying again']);
        end
    end
    
    if (attempt == max_attempts)
        add = false;
        disp('Failed to insert fibre, try increasing maximum number of iterations');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xc, tau] = fill_couette(o, rmin, rmax, nv, prams, seed)
   
    rng(seed);

    xc = zeros(2,nv);
    tau =zeros(1,nv);
    
    for i = 1:nv
       
        [Xnew, xc(:,i),tau(i), add] = o.add_fibre(rmin, rmax, prams,i); 
       
       if add
           o.X(:,end+1) = Xnew;
           
            oc = curve;
            o.center(:,end+1) = xc(:,i);
            o.nv = o.nv+1;
           [o.sa,o.xt,o.cur] = oc.diffProp(o.X);
           [~,o.length] = oc.geomProp(o.X);
       end

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = getXY(o)
    X = o.X;
end % getXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NearSelf,NearOther] = getZone(bd1,bd2,relate)
% [NearSelf,NearOther] = getZone(bd1,bd2,relate) constructs
% each boundary, index of the closest point, nearest point on a local
% interapolant, and argument of that nearest point.  bd1 contains
% the source points (which are also target points) and bd2
% contains additional target points.  The
% values of relate corresond to
% relate == 1  => only require NearSelf  (ie. bd1 to bd1)
% relate == 2  => only requirpe NearOther (ie. bd1 to bd2)
% relate == 3  => require both NearSelf and NearOther

if ~isempty(bd1.X)
    NearSelf = [];
    NearOther = [];
    
    N1 = bd1.N; % number of source/target points
    nv1 = bd1.nv; % number of source/target boundaries
    X1 = bd1.X; % source and target points
    oc = curve;
    [xsou,ysou] = oc.getXY(X1);
    % separate targets into x and y coordinates
    
    h = max(bd1.length)/N1;
    % smallest arclength over all boundaries
    ptsperbox = 10;
    % Estimate for number of points per box.  This simply sets the
    % number of uniformly refined boxes we take.  Estimate is not very
    % accurate.  What ptsperbox represents is the total number of points
    % that could be put in each two-dimensional bin where no two are
    % less than distance h from one another.  However, our points live
    % on curves and thus will not fill up an entire bin
    
    H = sqrt(ptsperbox)*h;
    xmin = min(min(xsou));
    xmax = max(max(xsou));
    xmin = xmin - H;
    xmax = xmax + H;
    ymin = min(min(ysou));
    ymax = max(max(ysou));
    ymin = ymin - H;
    ymax = ymax + H;
    % Add a buffer around the points so that it is easier to
    % work with bd2
    
    Nx = ceil((xmax - xmin)/H);
    Ny = ceil((ymax - ymin)/H);
    % Find bounds for box that contains all points and add a buffer
    % so that all points are guaranteed to be in the box
    
    Nbins = Nx * Ny; % Total number of bins
    
    ii = ceil((xsou - xmin)/H);
    jj = ceil((ysou - ymin)/H);
    % Index in x and y direction of the box containing each point
    bin = (jj-1)*Nx + ii;
    % Find bin of each point using lexiographic ordering (x then y)
    
    
    %figure(2);
    %clf; hold on
    %plot(xsou,ysou,'k.')
    %axis equal
    %axis([xmin xmin+Nx*H ymin ymin+Ny*H])
    %set(gca,'xtick',linspace(xmin,xmin+Nx*H,Nx+1))
    %set(gca,'ytick',linspace(ymin,ymin+Ny*H,Ny+1))
    %set(gca,'xticklabel',[]);
    %set(gca,'yticklabel',[]);
    %grid on
    %figure(1)
    %pause
    % DEBUG: This does a simple plot of the points with a grid that
    % aligns with the boundary of the boxes
    
    fpt = zeros(Nbins,nv1);
    lpt = zeros(Nbins,nv1);
    % allocate space for storing first and last points
    [binsort,permute] = sort(bin);
    % build permute.  Need binsort to find first and last points
    % in each box
    
    for k = 1:nv1 % Loop over boundaries
        for j = 1:N1 % Loop over bins
            ibox = binsort(j,k);
            if (fpt(ibox,k) == 0)
                fpt(ibox,k) = j;
                lpt(ibox,k) = 1;
            else
                lpt(ibox,k) = lpt(ibox,k) + 1;
            end
        end
        lpt(:,k) = fpt(:,k) + lpt(:,k) - 1;
    end
    % Construct first and last point in each box corresponding
    % to each boundary.  The order is based on permute.  For example,
    % permute(fpt(ibox,k)),...,permute(lpt(ibox,k)) is the set of
    % all points from boundary k contained in box ibox
    
    
    neigh = zeros(Nbins,9);
    
    %Do corners first
    neigh(1,1:4) = [1 2 Nx+1 Nx+2];
    % bottom left corner
    neigh(Nx,1:4) = [Nx Nx-1 2*Nx 2*Nx-1];
    % bottom right corner
    neigh(Nbins-Nx+1,1:4) = [Nbins-Nx+1 Nbins-Nx+2 ...
        Nbins-2*Nx+1 Nbins-2*Nx+2];
    % top left corner
    neigh(Nbins,1:4) = [Nbins Nbins-1 Nbins-Nx Nbins-Nx-1];
    % top right corner
    
    for j = 2:Nx-1
        neigh(j,1:6) = j + [-1 0 1 Nx-1 Nx Nx+1];
    end
    % neighbors of bottom row
    
    for j = Nbins-Nx+2:Nbins-1
        neigh(j,1:6) = j + [-1 0 1 -Nx-1 -Nx -Nx+1];
    end
    % neighbors of top row
    
    for j=Nx+1:Nx:Nbins-2*Nx+1
        neigh(j,1:6) = j + [-Nx -Nx+1 0 1 Nx Nx+1];
    end
    % neighbors of left column
    
    for j=2*Nx:Nx:Nbins-Nx
        neigh(j,1:6) = j + [-Nx-1 -Nx -1 0 Nx-1 Nx];
    end
    % neighbors of right column
    
    
    J = (Nx + 1:Nbins - Nx);
    J = J(mod(J-1,Nx)~=0);
    J = J(mod(J,Nx)~=0);
    % J is the index of boxes that are not on the boundary
    for j=J
        neigh(j,:) = j + [-Nx-1 -Nx -Nx+1 -1 0 1 Nx-1 Nx Nx+1];
    end
    % neighbors of interior points
    % TREE STRUCTURE IS COMPLETE
    
    
    if (relate == 1 || relate == 3)
        for k = 1:nv1
            distSS{k} = spalloc(N1,nv1,0);
            % dist(n,k,j) is the distance of point n on boundary k to boundary j
            zoneSS{k} = spalloc(N1,nv1,0);
            % near or far zone
            nearestSS{k} = spalloc(2*N1,nv1,0);
            % nearest point using local interpolant
            icpSS{k} = spalloc(N1,nv1,0);
            % index of closest discretization point
            argnearSS{k} = spalloc(N1,nv1,0);
            % argument in [0,1] of local interpolant
            nearFibersSS{k} = {};
        end
        % New way of representing near-singular integration structure so that
        % we can use sparse matricies.
        
        
        % begin classifying points where we are considering
        % boundary to boundary relationships
        for k = 1:nv1
            boxes = unique(bin(:,k));
            % Find all boxes containing points of boundary k
            boxes = neigh(boxes,:);
            % Look at all neighbors of boxes containing boundary k
            boxes = unique(boxes(:));
            % Remove repetition
            boxes = boxes(boxes~=0);
            % Delete non-existent boxes that came up because of neigh
            
            K = [(1:k-1) (k+1:nv1)];
            for k2 = K
                istart = fpt(boxes,k2);
                iend = lpt(boxes,k2);
                istart = istart(istart ~= 0);
                iend = iend(iend ~= -1);
                % Find index of all points in neighboring boxes of boundary k
                % that are in boundary k2
                
                neighpts = zeros(sum(iend-istart+1),1);
                % Allocate space to assign possible near points
                is = 1;
                for j=1:numel(istart)
                    ie = is + iend(j) - istart(j);
                    neighpts(is:ie) = permute(istart(j):iend(j),k2);
                    is = ie + 1;
                end
                % neighpts contains all points on boundary k2 that are in
                % neighboring boxes to boundary k
                
                neighpts = sort(neighpts);
                % sorting should help speedup as we won't be jumping around
                % through different boxes
                
                n = 0;
                for i=1:numel(neighpts)
                    ipt = neighpts(i);
                    ibox = bin(ipt,k2);
                    % box containing ipt on boundary k2
                    if (ibox ~= n)
                        n = ibox;
                        % Check if we need to move to a new box
                        neighbors = neigh(ibox,:);
                        % neighbors of this box
                        neighbors = neighbors(neighbors~=0);
                        % Remove non-existent neighbors
                        istart = fpt(neighbors,k);
                        iend = lpt(neighbors,k);
                        istart = istart(istart ~= 0);
                        iend = iend(iend ~= -1);
                        % Find points on boundary k in neighboring boxes
                        neighpts2 = zeros(sum(iend-istart+1),1);
                        is = 1;
                        for j=1:numel(istart)
                            ie = is + iend(j) - istart(j);
                            neighpts2(is:ie) = permute(istart(j):iend(j),k);
                            is = ie + 1;
                        end
                        % neighpts2 contains all points on boundary k that
                        % are in neighboring box of ibox
                    end % decide if we need to switch boxes
                    
                    [d0,d0loc] = min((xsou(ipt,k2) - xsou(:,k)).^2 + ...
                        (ysou(ipt,k2) - ysou(:,k)).^2);
                    % Find minimum distance between ipt on boundary k2 to
                    % possible closest points on boundary k
                    d0 = sqrt(d0);
                    % Save on not taking the square root on a vector but instead
                    % on a single real number
                    
                    icpSS{k}(ipt,k2) = d0loc;
                    if (d0 < 2*h);
                        [distSS{k}(ipt,k2),nearestx,nearesty,argnearSS{k}(ipt,k2)] = ...
                            bd1.closestPnt([xsou;ysou],xsou(ipt,k2),...
                            ysou(ipt,k2),k,icpSS{k}(ipt,k2));
                        nearestSS{k}(ipt,k2) = nearestx;
                        nearestSS{k}(ipt+N1,k2) = nearesty;
                        
                        nearFibersSS{k}= [nearFibersSS{k},k2];
                        % Find closest point along a local interpolant using
                        % Newton's method.
                        
                        if (distSS{k}(ipt,k2) < h)
                            zoneSS{k}(ipt,k2) = 1;
                        end
                        % Point ipt of boundary k2 is in the near zone of
                        % boundary k
                    end
                    
                    
                end % ipt
                
            end % k2
            
            nftmp = nearFibersSS{k};
            nearFibersSS{k} = unique([nftmp{:}]);
        end % k
        
        NearSelf.dist = distSS;
        NearSelf.zone = zoneSS;
        NearSelf.nearest = nearestSS;
        NearSelf.icp = icpSS;
        NearSelf.argnear = argnearSS;
        NearSelf.nearFibers = nearFibersSS;
        % Store everything in the structure NearSelf.  This way it is
        % much cleaner to pass everything around
        
    end % relate == 1 || relate == 3
    
    
    % Bin target points with respect to the source points
    if (relate == 2 || relate == 3)
        N2 = bd2.N; % number of additional targets
        nv2 = bd2.nv; % number of additional boundaries
        X2 = bd2.X; % additional target points
        [xtar,ytar] = oc.getXY(X2);
        % separate additional target points into x and y coordinates
        %  figure(2); clf
        %  plot(xtar,ytar,'r.')
        %  hold on
        %  plot(xsou,ysou,'k.')
        %  axis equal
        % DEBUG: FOR SEEING TARGET AND SOURCE POINTS IN THE TREE STRUCTURE
        % WHICH CAN BE PLOTTED ABOVE
        
        for k = 1:nv1
            distST{k} = spalloc(N1,nv2,0);
            % dist(n,k,j) is the distance of point n on boundary k to
            zoneST{k} = spalloc(N1,nv2,0);
            % near or far zone
            nearestST{k} = spalloc(2*N1,nv2,0);
            % nearest point using local interpolant
            icpST{k} = spalloc(N1,nv2,0);
            % index of closest discretization point
            argnearST{k} = spalloc(N1,nv2,0);
            % argument in [0,1] of local interpolant
        end
        % New way of representing near-singular integration structure so that
        % we can use sparse matricies.
        
        itar = ceil((xtar - xmin)/H);
        jtar = ceil((ytar - ymin)/H);
        [indx,indy] = find((itar >= 1) & (itar <= Nx) & ...
            (jtar >= 1) & (jtar <= Ny));
        % Only have to consider xx(ind),yy(ind) since all other points
        % are not contained in the box [xmin xmax] x [ymin ymax]
        
        for k = 1:nv1 % loop over sources
            for nind = 1:numel(indx)
                % loop over points that are not outside the box that surrounds
                % all target points with a sufficiently large buffer
                ii = indx(nind);
                jj = indy(nind);
                binTar = (jtar(ii,jj) - 1)*Nx + itar(ii,jj);
                boxesTar = neigh(binTar,:);
                boxesTar = boxesTar(boxesTar~=0);
                istart = fpt(boxesTar,k);
                iend = lpt(boxesTar,k);
                istart = istart(istart ~= 0);
                iend = iend(iend ~= -1);
                
                
                neighpts = zeros(sum(iend-istart+1),1);
                % Allocate space to assign possible near points
                if numel(neighpts) > 0
                    % it is possible of the neighboring boxes to contain
                    % no points.
                    is = 1;
                    for j = 1:numel(istart)
                        ie = is + iend(j) - istart(j);
                        neighpts(is:ie) = permute(istart(j):iend(j),k);
                        is = ie + 1;
                    end
                    % Set of potentially nearest points to
                    % (xtar(jj),ytar(jj))
                    
                    [d0,d0loc] = min((xtar(ii,jj) - xsou(neighpts,k)).^2 + ...
                        (ytar(ii,jj) - ysou(neighpts,k)).^2);
                    % find closest point and distance between (xtar(jj),ytar(jj))
                    % and boundary k.  Only need to look at points in neighboring
                    % boxes
                    d0 = sqrt(d0);
                    icpST{k}(ii,jj) = neighpts(d0loc);
                    
                    if d0 < 2*h
                        [distST{k}(ii,jj),nearestx,nearesty,argnearST{k}(ii,jj)] = ...
                            bd1.closestPnt([xsou;ysou],xtar(ii,jj),...
                            ytar(ii,jj),k,icpST{k}(ii,jj));
                        nearestST{k}(ii,jj) = nearestx;
                        nearestST{k}(ii+N2,jj) = nearesty;
                        
                        %          figure(2);
                        %          plot(xtar(ii,jj),ytar(ii,jj),'g.')
                        %          plot(nearestx,nearesty,'bo')
                        %          figure(1);
                        %          pause
                        % DEBUG: CHECK THAT NEWTON'S METHOD HAS DONE A GOOD JOB
                        % CONVERGING TO THE NEAREST POINT
                        % compute distance and nearest point between
                        % (xtar(ii,jj),ytar(ii,jj)) and boundary k
                        if distST{k}(ii,jj) < h
                            zoneST{k}(ii,jj) = 1;
                            % (xtar(ii,jj),ytar(ii,jj)) is in the near zone of boundary k
                        end
                    end % d0 < 2*h
                end % numel(neighpts) > 0
                
            end % ii and jj
            
        end % k
        
        NearOther.dist = distST;
        NearOther.zone = zoneST;
        NearOther.nearest = nearestST;
        NearOther.icp = icpST;
        NearOther.argnear = argnearST;
        NearOther.nearFibers = [];
        % store near-singluar integration requirements in structure NearOther
        
    end % relate == 2 || relate == 3
else
    NearSelf  = [];
    NearOther = [];
    
    NearSelf.dist = [];
    NearSelf.zone = [];
    NearSelf.nearest = [];
    NearSelf.icp = [];
    NearSelf.argnear = [];
    NearSelf.nearFibers = [];
    
    NearOther.dist = [];
    NearOther.zone = [];
    NearOther.nearest = [];
    NearOther.icp = [];
    NearOther.argnear = [];
    NearOther.nearFibers = [];
end
end % getZone


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist,nearestx,nearesty,theta] = closestPnt(o,X,xTar,ytar,k,icp)
% [dist,nearestx,nearesty,theta] = closestPnt(X,xTar,ytar,k,icp)
% computes the closest point on boundary k to (xTar,ytar)
% using a Lagrange interpolant.  icp is the index of the closest
% point on the discrete mesh which is used as an initial guess

N = size(X,1)/2; % Number of points per boundary 
A = poten.lagrangeInterp;
interpOrder = size(A,1);
% need interpolation matrix and its size

p = ceil((interpOrder+1)/2);
% Accommodate for either an even or odd number of interpolation
% points
pn = mod((icp-p+1:icp-p+interpOrder)' - 1,N) + 1;
% band of points around icp.  The -1,+1 combination sets index
% 0 to N as required by the code

px = A*X(pn,k); % polynomial interpolant of x-coordinate
py = A*X(pn+N,k); % polynomial interpolant of y-coordinate
Dpx = px(1:end-1).*(interpOrder-1:-1:1)';
Dpy = py(1:end-1).*(interpOrder-1:-1:1)';
D2px = Dpx(1:end-1).*(interpOrder-2:-1:1)';
D2py = Dpy(1:end-1).*(interpOrder-2:-1:1)';
% To do Newton's method, need two derivatives

theta = 1/2;
% midpoint is a good initial guess
for newton = 1:2
  zx = filter(1,[1 -theta],px);
  zx = zx(end);
  zy = filter(1,[1 -theta],py);
  zy = zy(end);
  Dzx = filter(1,[1 -theta],Dpx);
  Dzx = Dzx(end);
  Dzy = filter(1,[1 -theta],Dpy);
  Dzy = Dzy(end);
  D2zx = filter(1,[1 -theta],D2px);
  D2zx = D2zx(end);
  D2zy = filter(1,[1 -theta],D2py);
  D2zy = D2zy(end);
  % Using filter is the same as polyval, but it is much
  % faster when only requiring a single polyval such as here.

  newtonNum = (zx-xTar)*Dzx + (zy-ytar)*Dzy;
  % numerator of Newton's method
  newtonDen = (zx-xTar)*D2zx + (zy-ytar)*D2zy + ...
      Dzx^2 + Dzy^2;
  % denominator of Newton's method
  theta = theta - newtonNum/newtonDen;
  % one step of Newton's method
end
% Do a few (no more than 3) Newton iterations

nearestx = filter(1,[1,-theta],px);
nearestx = nearestx(end);
nearesty = filter(1,[1,-theta],py);
nearesty = nearesty(end);
dist = sqrt((nearestx - xTar)^2 + (nearesty - ytar)^2);
% Compute nearest point and its distance from the target point

end % closestPnt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function icollision = collision(geom,near,ifmm,inear,op, om)

if (om.profile)
    tCollision = tic;
end

f = [ones(geom.N,geom.nv);zeros(geom.N,geom.nv)];
% density function.  Solving a scalar-valued layer-potential, so set the
% second component of density function to 0

if (ifmm)
    kernel = @op.exactLaplaceDLfmm;
else
    kernel = @op.exactLaplaceDL;
end

D = @(X) zeros(2*size(X,1),size(X,2));
% know that the limiting value of the DLP of the density function is
% always zero, so don't have to build the DLP matrix for
% self-interactions
if inear
  Fdlp = op.nearSingInt(geom,f,D,[],near,kernel,kernel,geom,true,false);
else
  Fdlp = kernel(geom,f,D);
end

Fdlp = Fdlp(1:geom.N,:);
% take only first component of the DLP
% plot(Fdlp);

buffer = 1e-4;

icollision = any(abs(Fdlp(:)) > buffer);

if (om.profile)
    om.writeMessage(['Collision detection completed in ', num2str(toc(tCollision)), ' seconds']);
end

end
% collision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function press = pressure(geom,f,stokeslets,pressTar,sEqualsT)
% press = pressure(vesicle,f,RS,pressTar,fmm,LP) computes the pressure
% due to vesicle at the locations pressTar with near-singular
% integration with traction jump f. Note that vesicle may correspond to 
% solid walls rather than vesicles

N = geom.N; % points per vesicle
nv = geom.nv; % number of vesicles
Nup = N*ceil(sqrt(N));
op = poten(N);
Ntar = pressTar.N;
% number of target points where we want to compute the pressure
% and stress

if sEqualsT
    [NearStruct,~] = geom.getZone(pressTar,1);
else
    [~,NearStruct] = geom.getZone(pressTar,2);
end
% build near-singular integration structure for vesicle
% to pressure target points
%InOutFlag = find(vesicle.sortPts(pressTar.X,fmm,NearV2T) == 1);
%figure(1); clf; hold on
%plot(vesicle.X(1:end/2,:),vesicle.X(end/2+1:end,:),'r');
%oc = curve;
%[xx,yy] = oc.getXY(pressTar.X);
%plot(xx,yy,'b.');
%plot(xx(InOutFlag),yy(InOutFlag),'k.');
% DEBUG: MAKE SURE POINTS INSIDE AND OUTSIDE ARE LABELLED CORRECTLY

oc = curve;
PdiagIn = op.pressDLmatrix(geom);
PdiagOut = PdiagIn;
Deriv = fft1.D1(N);
% spectral differentiation matrix
for k = 1:nv
    [tanx,tany] = oc.getXY(geom.xt(:,k));
    % tangent vector
    sa = geom.sa(:,k);
    % Jacobian
    
    jump = -[diag(tanx./sa) diag(tany./sa)] * ...
        [Deriv zeros(N); zeros(N) Deriv];
    PdiagIn(:,:,k) = PdiagIn(:,:,k) - jump;
    PdiagOut(:,:,k) = PdiagOut(:,:,k) + jump;
    % add in the jump term.  Jump term is negative because target
    % points are interior to the solid walls
end
% Jump term that comes from pressure of double-layer potential
% is the dot product of the tangential derivative of the traction
% jump with the tangent.  Assuming that all the pressure target
% points are inside the physical domain so there is no use to
% consider the limiting value from the other side
Pdiag = @(X) op.exactPressureDLdiag(geom,[PdiagOut;PdiagIn],X);
% nearSingInt assumes that input and output are vector-valued Take
% adavntage of the fact that nearSingInt only works with
% vector-valued densities by passing in two different jumps---one for
% the limit from the inside and one from the outside.

kernel = @op.exactPressDL;
kernelDirect = @op.exactPressDL;

Xup = [interpft(geom.X(1:N,:),Nup); interpft(geom.X(N+1:2*N,:),Nup)];

geomUp = capsules([],Xup);
Dup = op.stokesDLmatrix(geomUp);

press = op.nearSingInt(geom,f,Pdiag,Dup,NearStruct,kernel,kernelDirect,...
                    pressTar,sEqualsT,false);

% compute the pressure of the single- or double-layer 
% potential using near-singular integration.  First row 
% corresponds to using the limiting value for the pressure 
% exterior of the vesicle and the second row corresponds to 
% using the limiting value for the pressure interior of the 
% vesicle
% press = pressT(1:Ntar);
press = press(1:Ntar,:);
%press(InOutFlag) = pressT(Ntar+InOutFlag);
% At the interior points, take the second row

if sEqualsT
    % add self contribution
    for k = 1:geom.nv        
        press(:,k) = press(:,k) + PdiagIn(:,:,k)*f(:,k);
    end
end


nwalls = length(stokeslets)/2; % number of inner walls

if nwalls > 0
  for k = 1:nwalls
    cx = geom.center(1,k+1);
    cy = geom.center(2,k+1);
    % center of interior solid wall k
    xi1 = stokeslets(2*(k-1)+1);
    xi2 = stokeslets(2*(k-1)+2);
    % first and second components of the stokeslet on 
    % interior solid wall k
    [x,y] = oc.getXY(pressTar.X);
    rx = x - cx;
    ry = y - cy;
    rho2 = rx.^2 + ry.^2;
    press = press +  2./(rx.^2 + ry.^2).*(rx*xi1 + ry*xi2);
    % add in pressure of stokeslet term
    % rotlet has a vanishing pressure gradient
  end
end

end % pressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stress1,stress2] = stressTensor(vesicle,f,RS,stressTar,sEqualsT)
% [stress1 stress2] = stressTensor(vesicle,f,RS,stressTar,fmm)
% computes the stress tensor due to vesicle at the locations stressTar
% with or without near-singular integration with traction jump f.
% Returns in two components.  stress1 is the stress tensor applied to
% [1;0] and stress2 is the stress tensor applied to [0;1].
% Note that vesicle may correspond to solid walls rather
% than vesicles.  If doing solid walls, there is a stress due to the
% rotlet and stokeslet terms that needs to be added in

N = vesicle.N; % points per vesicle
nv = vesicle.nv; % number of vesicles
op = poten(N);
Ntar = stressTar.N*stressTar.nv;
% number of points where we want to evaluate the stress

tangent = vesicle.xt;
oc = curve;
[tanx,tany] = oc.getXY(tangent);
nx = tany;
ny = -tanx;
% decompose tangent and normal vectors into their x and y
% coordinates

if sEqualsT
    [NearStruct,~] = vesicle.getZone(stressTar,1);
else
    [~,NearStruct] = vesicle.getZone(stressTar,2);
end
% build near-singular integration structure for vesicle
% to pressure target points

    
[S1diagIn,S2diagIn] = op.stressDLmatrix(vesicle);
% use odd-even integration to evaluate stress due to
% each vesicle independent of the others on its boundary
% S1diagIn is the matrix that takes the traction jump and returns
% the stress tensor applied to [1;0] evaluated on the boundary of
% the vesicle
% S2diagIn is the matrix that takes the traction jump and returns
% the stress tensor applied to [0;1] evaluated on the boundary of
% the vesicle

S1diagOut = S1diagIn;
S2diagOut = S2diagIn;
% Use a second copy for the jump from the outside
Deriv = fft1.D1(N);
% Fourier differentiation matrix

for k = 1:nv
    [tanx,tany] = oc.getXY(vesicle.xt(:,k));
    % tangent vector
    sa = vesicle.sa(:,k);
    % Jacobian

    jump = diag(1./sa.*(1+tanx.^2-tany.^2).*tanx)*Deriv;
    S1diagIn(1:N,1:N,k) = S1diagIn(1:N,1:N,k) - jump;
    S1diagOut(1:N,1:N,k) = S1diagOut(1:N,1:N,k) + jump;
    % add in (1,1) jump component to stress applied to [1;0]

    jump = diag(1./sa.*(1+tanx.^2-tany.^2).*tany)*Deriv;
    S1diagIn(1:N,N+1:2*N,k) = S1diagIn(1:N,N+1:2*N,k) - jump;
    S1diagOut(1:N,N+1:2*N,k) = S1diagOut(1:N,N+1:2*N,k) + jump;
    % add in (1,2) jump component to stress applied to [1;0]

    jump = diag(1./sa.*(2*tanx.*tany).*tanx)*Deriv;
    S1diagIn(N+1:2*N,1:N,k) = S1diagIn(N+1:2*N,1:N,k) - jump;
    S1diagOut(N+1:2*N,1:N,k) = S1diagOut(N+1:2*N,1:N,k) + jump;
    % add in (2,1) jump component to stress applied to [1;0]

    jump = diag(1./sa.*(2*tanx.*tany).*tany)*Deriv;
    S1diagIn(N+1:2*N,N+1:2*N,k) = ...
        S1diagIn(N+1:2*N,N+1:2*N,k) - jump;
    S1diagOut(N+1:2*N,N+1:2*N,k) = ...
        S1diagOut(N+1:2*N,N+1:2*N,k) + jump;
    % add in (2,2) jump component to stress applied to [1;0]

    jump = diag(1./sa.*(2*tanx.*tany).*tanx)*Deriv;
    S2diagIn(1:N,1:N,k) = S2diagIn(1:N,1:N,k) - jump;
    S2diagOut(1:N,1:N,k) = S2diagOut(1:N,1:N,k) + jump;
    % add in (1,1) jump component to stress applied to [0;1]

    jump = diag(1./sa.*(2*tanx.*tany).*tany)*Deriv;
    S2diagIn(1:N,N+1:2*N,k) = S2diagIn(1:N,N+1:2*N,k) - jump;
    S2diagOut(1:N,N+1:2*N,k) = S2diagOut(1:N,N+1:2*N,k) + jump;
    % add in (1,2) jump component to stress applied to [0;1]

    jump = diag(1./sa.*(1-tanx.^2+tany.^2).*tanx)*Deriv;
    S2diagIn(N+1:2*N,1:N,k) = S2diagIn(N+1:2*N,1:N,k) - jump;
    S2diagOut(N+1:2*N,1:N,k) = S2diagOut(N+1:2*N,1:N,k) + jump;
    % add in (2,1) jump component to stress applied to [0;1]

    jump = diag(1./sa.*(1-tanx.^2+tany.^2).*tany)*Deriv;
    S2diagIn(N+1:2*N,N+1:2*N,k) = ...
        S2diagIn(N+1:2*N,N+1:2*N,k) - jump;
    S2diagOut(N+1:2*N,N+1:2*N,k) = ...
        S2diagOut(N+1:2*N,N+1:2*N,k) + jump;
    % add in (2,2) jump component to stress applied to [0;1]
end
% Jump term that comes from stress of double-layer potential
S1diagOutFn = @(X) op.exactStressDLdiag(vesicle,S1diagOut,X);
S2diagOutFn = @(X) op.exactStressDLdiag(vesicle,S2diagOut,X);
% nearSingInt assumes that input and output are vector-valued

kernel1 = @op.exactStressDL1;
kernel2 = @op.exactStressDL2;

% Have built all single- or double-layer potentials on boundary so that
% we can compute diagonal value to use near-singular integration

% CREATE UPSAMPLED MATRICES
Xsou = vesicle.X; 
Nup = N*ceil(sqrt(N));

Xup = [interpft(Xsou(1:N,:),Nup);...
       interpft(Xsou(N+1:2*N,:),Nup)];

geomUp = capsules([],Xup);
[Dup1, Dup2] = op.stressDLmatrix(geomUp);

stress1 = op.nearSingInt(vesicle,f,S1diagOutFn,Dup1,NearStruct,kernel1,kernel1,...
    stressTar,sEqualsT,false);
  % correct stress at exterior points

% Use near-singular integration to compute first component of
% the stress due to the double-layer potential

stress2 = op.nearSingInt(vesicle,f,S2diagOutFn,Dup2,NearStruct,kernel2,kernel2,...
    stressTar,sEqualsT,false);
% correct stress at exterior points

if sEqualsT
    % add self contribution
    for k = 1:vesicle.nv        
        stress1(:,k) = stress1(:,k) + S1diagOut(:,:,k)*f(:,k);
        stress2(:,k) = stress2(:,k) + S2diagOut(:,:,k)*f(:,k);
    end
end
% Use near-singular integration to compute second component of
% the stress due to the double-layer potential

if ~isempty(RS)
  for k = 2:vesicle.nv
    cx = vesicle.center(1,k);
    cy = vesicle.center(2,k);
    % center of interior solid wall k
    xi = RS(3*(k-2) + 3);
    % rotlet on interior solid wall k
    xi1 = RS(3*(k-2) + 1);
    xi2 = RS(3*(k-2) + 2);
    % first and second components of the stokeslet on 
    % interior solid wall k 
    [x,y] = oc.getXY(stressTar.X);
    rx = x - cx;
    ry = y - cy;
    rho4 = (rx.^2 + ry.^2).^2;
    
    stress1(1:stressTar.N,:) = stress1(1:stressTar.N,:) - 4*xi*rx.*ry./rho4;
    % add in (1,1) component of stress due to Rotlet
    stress1(1+stressTar.N:end,:) = stress1(1+stressTar.N:end,:) + ...
                    2*xi*(rx.^2 - ry.^2)./rho4;
    % add in (1,2) component of stress due to Rotlet
    stress2(1:stressTar.N,:) = stress2(1:stressTar.N,:) + 2*xi*(rx.^2 - ry.^2)./rho4;
    % add in (2,1) component of stress due to Rotlet
    stress2(1+stressTar.N:end,:) = stress2(1+stressTar.N:end,:) + 4*xi*rx.*ry./rho4;
    % add in (2,2) component of stress due to Rotlet


    stress1(1:stressTar.N,:) = stress1(1:stressTar.N,:) + ...
                2*(rx*xi1 + ry*xi2).*(-rx.^2 + ry.^2)./rho4;
    stress1(1:stressTar.N,:) = stress1(1:stressTar.N,:) - ...
                2*(rx*xi1 + ry*xi2)./(rx.^2+ry.^2);
    % add in (1,1) component of stress due to Stokeslet
    % need to add in -1*pressure term as well
    stress1(1+stressTar.N:end,:) = stress1(1+stressTar.N:end,:) - ...
                4*rx.*ry.*(rx*xi1 + ry*xi2)./rho4;
    % add in (1,2) component of stress due to Stokeslet
    stress2(1:stressTar.N,:) = stress2(1:stressTar.N,:) - ...
                4*rx.*ry.*(rx*xi1 + ry*xi2)./rho4;
    % add in (1,2) component of stress due to Stokeslet
    stress2(1+stressTar.N:end,:) = stress2(1+stressTar.N:end,:) + ...
                2*(rx*xi1 + ry*xi2).*(rx.^2 - ry.^2)./rho4;
    stress2(1+stressTar.N:end,:) = stress2(1+stressTar.N:end,:) - ...
                2*(rx*xi1 + ry*xi2)./(rx.^2+ry.^2);
    % add in (2,2) component of stress due to Stokeslet
    % need to add in -1*pressure term as well
  end
end
% Add in contributions due to Rotlets and Stokeslets

end % stressTensor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InOut = sortPts(vesicle,Xtar,fmm,NearV2T)
% InOut = sortPts(vesicle,Xtar,fmm,nearV2T) determines if the set of
% points in Xtar are inside or outside of a vesicle configuration

density = [ones(vesicle.N,vesicle.nv);...
           zeros(vesicle.N,vesicle.nv)];
% density function that is used to check if a point is inside or
% outside

tracers = capsules([], Xtar);
if nargin == 3 
  [~,NearV2T] = vesicle.getZone(tracers,2);
end

op = poten(vesicle.N);
if ~fmm
  kernel = @op.exactLaplaceDL;
else
  kernel = @op.exactLaplaceDLfmm;
end
% kernel for Laplace's double layer potential

DLP = @(X) zeros(2*size(X,1),size(X,2));
% can cheat here because we know that the double-layer potential
% applied to the constant density function will always return zero

kernelDirect = kernel;
InOut = op.nearSingInt(vesicle,density,DLP,[],NearV2T,kernel,kernelDirect,tracers,false,false);

InOut = InOut(1:end/2);
% only care about the first half as we are working with laplace as
% opposed to the vector-valued stokes layer potentials

thresh = 1e-7;
InOut(abs(InOut) > thresh) = 1;
InOut(abs(InOut) <= thresh) = 0;
% for points inside a vesicle, but close to its boundary, we will be
% interpolating the function [0;1;1;1;1...] close to the initial
% point.  Therefore, the value returned by the DLP will be some value
% between 0 and 1, but it can be quite close to 0.  So, we threshold

end % sortPts

end % methods



end %capsules



