classdef post
% Post-processing of data from output files

properties    
prams
options

times
xc
tau
U
omega
eta_p
eta_w
rotlets
stokeslets
wall_c

OUTPUTPATH_GIFS;

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = post(data_file)
 
load(data_file)

o.prams = prams;
o.options = options;
o.times = t;
o.xc = xc;
o.tau = tau;
o.U = U;
o.omega = omega;
o.eta_p = eta_p;
o.eta_w = eta_w;
o.rotlets = rot;
o.stokeslets = stokes;
%o.wall_c = wall_c;

o.OUTPUTPATH_GIFS = '../output/gifs/';

end %post : constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotFibres(o, iT)
    
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
    fill3(geom.X(1:end/2,:),geom.X(end/2+1:end,:),...
                        100*ones(o.prams.Np,o.prams.np), 'k');
    
end % post : plot_fibres

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotWalls(o,iT)
    
    oc = curve;
    xWalls = oc.createWalls(o.prams.Nw, o.options);
    
    hold on;
    for i = 1:o.prams.nw
        plot3([xWalls(1:end/2,i);xWalls(1,i)], ...
            [xWalls(end/2+1:end,i);xWalls(end/2+1,i)],...
            99*ones(o.prams.Nw+1, 1),'r', 'linewidth', 2);
    end

    ptsTrack = o.prams.tracker_fnc(o.times(iT));
    
    plot3(ptsTrack(:,1),ptsTrack(:,2), 100*ones(size(ptsTrack,1),1),...
                    'b.', 'MarkerSize', 20); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = animatedGif(o, gif_options)
    
h = figure();
set(h,'Units','normalized');
set(h,'Position',[0 0 1 1]);

[~,~,m] = size(o.xc);
if strcmp(gif_options.itmax, 'all')
    itmax = m;
else
    itmax = gif_options.itmax;
end

% get axis limits
if strcmp(gif_options.xmin, 'auto:all')
    xmin = min(min(o.xc(1,:,1:itmax))) - o.prams.lengths;
else
    if ~strcmp(gif_options.xmin, 'auto:frame')
        xmin = gif_options.xmin;
    end
end

if strcmp(gif_options.xmax, 'auto:all')
    xmax = max(max(o.xc(1,:,1:itmax))) + o.prams.lengths;
else
    if ~strcmp(gif_options.xmax, 'auto:frame')
        xmax = gif_options.xmax;
    end
end

if strcmp(gif_options.ymin, 'auto:all')
    ymin = min(min(o.xc(2,:,1:itmax))) - o.prams.lengths;
else
    if ~strcmp(gif_options.ymin, 'auto:frame')
        ymin = gif_options.ymin;
    end
end

if strcmp(gif_options.ymax, 'auto:all')
    ymax = max(max(o.xc(2,:,1:itmax))) + o.prams.lengths;
else
    if ~strcmp(gif_options.ymax, 'auto:frame')
        ymax = gif_options.ymax;
    end
end
   
% loop over all time steps
for i = 1:gif_options.stride:itmax-1
    
    clf;
    axis equal
    
    hold on
    
    if ~gif_options.axis 
        axis off
    end
    
    if strcmp(gif_options.xmin, 'auto:frame')
        xmin = min(min(o.xc(1,:,i))) - max(o.prams.lengths);
    end
    
    if strcmp(gif_options.xmin, 'auto:frame')
        xmax = max(max(o.xc(1,:,i))) + max(o.prams.lengths);    
    end
    
    if strcmp(gif_options.xmin, 'auto:frame')
        ymin = min(min(o.xc(2,:,i))) - max(o.prams.lengths);
    end
    
    if strcmp(gif_options.xmin, 'auto:frame')
        ymax = max(max(o.xc(2,:,i))) + max(o.prams.lengths);
    end
      
    if ~isempty(gif_options.contour_field)
        couette_u = @(x,y) -(y.*(100-x.^2-y.^2))./(99*(x.^2+y.^2));
        if (i > 1)
            o.plotContourCircle(gif_options.contour_field, i, 5.01, 9.99, 0, 0, 100,  300, ...
                                false, cmin, cmax, couette_u);
        else
            [~,~,~,cmin,cmax] = o.plotContourCircle(gif_options.contour_field, i, 5.01, 9.99, ...
                                    0, 0, 100, 300, false, [], [], couette_u);
        end
        
        %o.plotContourBox('pressure', i, -2, 10, -3, 3, 300, 100, false, @(x,y) 0);
    end
    
    if gif_options.velocity_quiver
         o.velocityQuiverBox(i, xmin, xmax, ymin, ymax, 20, 20, 1, gif_options.bg_flow);
    end

    if o.options.confined
        o.plotWalls(i);
    end
    
    if o.prams.np > 0
        o.plotFibres(i);
    end
    
    axis equal
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
   
    title(sprintf('t = %6.3f', o.times(i)));
    
    drawnow;
    pause(0.1);
    frame = getframe(h);
    im = frame2im(frame);
    
    [imind, cm] = rgb2ind(im,256);
    
    % output as sequence of files
    
    if strcmp(gif_options.file_type, 'tikz')
        % output as tex files
        addpath('../src/matlab2tikz/src');
        
        if i == 1
            mkdir([o.OUTPUTPATH_GIFS, gif_options.file_name]);
        end
        
        cleanfigure();
        matlab2tikz([o.OUTPUTPATH_GIFS, gif_options.file_name, '/', ...
                 gif_options.file_name, '-', sprintf('%03d', i), '.tikz'],...
                 'height', '10cm', 'width', '12cm', 'standalone', true, 'floatFormat', '%.3f');
        
    else if strcmp(gif_options.file_type, 'gif')
            if i == 1;
                
                imwrite(imind, cm, ...
                    [o.OUTPUTPATH_GIFS,  gif_options.file_name, '.gif'], ...
                    'gif', 'Loopcount',inf, 'DelayTime', 0.1);
                
            else
                
                imwrite(imind, cm,...
                    [o.OUTPUTPATH_GIFS,  gif_options.file_name, '.gif'], ...
                    'gif','WriteMode','append', 'DelayTime',0.1);               
                
            end
        end
    end
end

end % post : aimated_gif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = evaluateDLP(~, geom, eta, X)
% evaluates the Stokes double layer potential at target points X, given a 
% geometry geom and a density function eta


geomTar = capsules([],X);

pot = poten(geom.N);
[~,NearStruct] = geom.getZone(geomTar,2);

D = pot.stokesDLmatrix(geom);
DLP = @(X) pot.exactStokesDLdiag(geom,D,X) - 1/2*X;
u = pot.nearSingInt(geom, eta, DLP,[],...
    NearStruct, @pot.exactStokesDLfmm, @pot.exactStokesDL, geomTar,false,false);

end % post : evaluateDLP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = evaluateRotStokes(o, iT, X)
    
xi = o.rotlets(:,iT);
lambda = o.stokeslets(:,:,iT);

oc = curve;
[xWalls,cen] = oc.createWalls(o.prams.Nw, o.options);
walls = capsules([],xWalls);
walls.center = cen;

[x,y] = oc.getXY(X);

u = zeros(2,length(x));

for i = 1:length(xi)
    r = [x - cen(i+1,1); y - cen(i+1,2)];
    %u = u + xi(i)*[r(2,:); -r(1,:)]; % add contributation from rotlets

    for j = 1:length(x)
       
        rho = norm(r(:,j));
        u(:,j) = u(:,j) + xi(i)*[r(2,j);-r(1,j)]/rho^2; % add contributation from rotlets

        ror = r(:,j)*r(:,j)';

        u(:,j) = u(:,j) + (ror/rho^2 - log(rho)*eye(2))*lambda(:,i);
    end
end
end % post : evaluateRotStokes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = evaluateVelocity(o, iT, X)

u = zeros(size(X));

if o.prams.np > 0 % add contribution from fibers
    
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
    u = u + o.evaluateDLP(geom, o.eta_p(:,:,iT), X);
end

if o.options.confined
    oc = curve;
    xWalls = oc.createWalls(o.prams.Nw, o.options);
    walls = capsules(o.prams,xWalls);
    u = u + o.evaluateDLP(walls, o.eta_w(:,:,iT), X);
    
    if o.prams.nw > 1
        u = u + o.evaluateRotStokes(iT, X);
    end
end

end % post : evaluateVelocity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pressure = evaluatePressure(o, iT, X)
% evaluates the stress tensor at time step iT at target points X, given a 
% boundary type, either 'fibers', or 'walls'

if o.prams.np > 0
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
else
    geom = capsules(o.prams, [], []);
end

geomTar = capsules([],X);

if geom.n > 0
    pressureF = geom.pressure(o.eta_p(:,:,iT),[],geomTar,false);
else
    pressureF = zeros(1,length(X));
end

if o.options.confined
    oc = curve;
    xWalls = oc.createWalls(o.prams.Nw, o.options);
    walls = capsules(o.prams,xWalls);
    
    if o.prams.nw > 1
        S = o.stokeslets(:,:,iT);
    else
        S = [];
    end
    pressureW = walls.pressure(o.eta_w(:,:,iT),S,geomTar,false);
else
    pressureW = zeros(1, length(X));
end

pressure = pressureF + pressureW;

end % post : evaluatePressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[stress1, stress2] = evaluateStress(o, iT, X)
% evaluates the stress tensor at time step iT at target points X, given a
% boundary type, either 'fibers', or 'walls'

if o.prams.np > 0
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
else
    geom = capsules(o.prams, [], []);
end

geomTar = capsules([],X);

RS = [];
if geom.n > 0
    [stress1,stress2] = geom.stressTensor(o.eta_p(:,:,iT),RS,geomTar,false);
else
    stress1 = zeros(2,geomTar.n);
    stress2 = zeros(2,geomTar.n);
end

if o.options.confined

    RS = zeros(3*(o.prams.w),1);
    for k = 1:o.prams.nw-1
        RS(3*(k-1)+1:3*(k-1)+2) = o.stokeslets(2*k-1:2*k);
        RS(3*(k-1)+3) = o.rotlets(k);
    end
    
    oc = curve;
    [xWalls,cen] = oc.createWalls(o.prams.Nw, o.options);
    walls = capsules([],xWalls);
    
    walls.center = cen;
    
    [stress1w,stress2w] = walls.stressTensor(o.eta_w(:,:,iT),RS,geomTar,false);
    
    stress1 = stress1 + stress1w;
    stress2 = stress2 + stress2w;
end

end % post : evaluateStress

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dissipation = evaluateDissipation(o, iT, X)

[stress1,stress2] = o.evaluateStress(iT, X);
pressure = o.evaluatePressure(iT,X);

dissipation = zeros(length(X),1);

for i = 1:length(dissipation)
    stress_visc = [stress1(:,i), stress2(:,i)] + eye(2)*pressure(i);
    dissipation(i) = 0.5*norm(stress_visc,'fro')^2;
end

end % post : evaluateDissipation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,U,V] = velocityQuiver(o, iT, X, Y, stride, bg_flow)

[M,N] = size(X);
    
oc = curve;
if o.prams.np > 0
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
else
    geom = capsules(o.prams, [], []);
end

W = o.evaluateVelocity(iT, [X(:)';Y(:)']);
U = W(1,:);
V = W(2,:);


%add in background flow
bg = bg_flow(X(:),Y(:))';
U = U - bg(1,:);
V = V - bg(2,:);

if o.prams.np > 0
    insideFibers = geom.sortPts([X(:);Y(:)],true);
    W(insideFibers==1) = NaN;
    X(insideFibers==1) = NaN;
    Y(insideFibers==1) = NaN;
end

if o.options.confined
    xWalls = oc.createWalls(o.prams.Nw, o.options);
    
    
    wallOuter = capsules(o.prams, xWalls(:,1));
    insideOuterWall = wallOuter.sortPts([X(:);Y(:)],true);
    
    % identify points outside outer wall
    W(insideOuterWall==0) = NaN;
    
    if o.prams.nw > 1
        wallsInner = capsules(o.prams, xWalls(:,2:end));
        
        insideInnerWalls = wallsInner.sortPts([X(:);Y(:)],true);
        % identify points inside inner walls
        W(insideInnerWalls==1) = NaN;
        X(insideInnerWalls==1) = NaN;
        Y(insideInnerWalls==1) = NaN;
    end
    
    X(insideOuterWall==0) = NaN;
    Y(insideOuterWall==0) = NaN;
end    
 
U = reshape(U,M,N);
V = reshape(V,M,N);


quiver(X, Y, U, V, stride);

end % post : velocityQuiver

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,U,V] = velocityQuiverBox(o, iT, xmin, xmax, ymin, ymax, ...
                        Nx, Ny, stride, bg_flow)
    
[X, Y] = meshgrid(linspace(xmin, xmax, Nx), linspace(ymin, ymax, Ny));
[X,Y,U,V] = o.velocityQuiver(iT, X, Y, stride, bg_flow);

end % post : plotContourBox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,W, cmin, cmax] = plotContour(o, variable, iT, X, Y, ...
                        replace_inside, cmin, cmax, exact_sol)

[M,N] = size(X);
    
oc = curve;
if o.prams.np > 0
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
else
    geom = capsules(o.prams, [], []);
end

switch variable

    case 'velocity_u'
     
        W = o.evaluateVelocity(iT, [X(:)';Y(:)']);
        W = W(1,:);
        
    case 'velocity_v'
        
        W = o.evaluateVelocity(iT,  [X(:)';Y(:)']);
        W = W(2,:);
        
    case 'velocity_u_err'
        
        W = o.evaluateVelocity(iT, [X(:)';Y(:)']);
        W = W(1,:);
        W_exact = exact_sol(X(:)',Y(:)');
        W = W - W_exact;
        
    case 'velocity_v_err'
        
        W = o.evaluateVelocity(iT,  [X(:)';Y(:)']);
        W = W(2,:);   
        W_exact = exact_sol(X(:)',Y(:)');
        W = W - W_exact;
        
    case 'pressure'
        W = o.evaluatePressure(iT,  [X(:)';Y(:)']);
        
    case 'pressure_err'
        
        W = o.evaluatePressure(iT,  [X(:)';Y(:)']);
        W_exact = exact_sol(X(:)',Y(:)');
        W = W - W_exact;
        
    case 'stress_xx'
        
        [W, ~] = o.evaluateStress(iT,  [X(:)';Y(:)']);
        W = W(1,:);
        
    case 'stress_xy'
        
        [W, ~] = o.evaluateStress(iT,  [X(:)';Y(:)']);
        W = W(2,:);
        
    case 'stress_yx'
        
        [~, W] = o.evaluateStress(iT,  [X(:)';Y(:)']);
        W = W(1,:);
        
    case 'stress_yy'
        
        [~, W] = o.evaluateStress(iT,  [X(:)';Y(:)']);
        W = W(2,:);
        
    case 'dissipation'
        
        W = o.evaluateDissipation(iT,  [X(:)';Y(:)']);
end


if replace_inside
    % indentify points inside fiber
    insideFibers = geom.sortPts([X(:);Y(:)],true);
    
    if o.options.confined
        xWalls = oc.createWalls(o.prams.Nw, o.options);
        
        
        wallOuter = capsules(o.prams, xWalls(:,1));
        insideOuterWall = wallOuter.sortPts([X(:);Y(:)],true);
        
        % identify points outside outer wall
        W(insideOuterWall==0) = NaN;
        
        if o.prams.nw > 1
            wallsInner = capsules(o.prams, xWalls(:,2:end));
            
            insideInnerWalls = wallsInner.sortPts([X(:);Y(:)],true);
            % identify points inside inner walls
            W(insideInnerWalls==1) = NaN;
            X(insideInnerWalls==1) = NaN;
            Y(insideInnerWalls==1) = NaN;
        end
        
        X(insideOuterWall==0) = NaN;
        Y(insideOuterWall==0) = NaN;
    end
    
    W(insideFibers==1) = NaN;
    X(insideFibers==1) = NaN;
    Y(insideFibers==1) = NaN;
end
W = reshape(W,M,N);

surf(X,Y,W);
view(2);
shading interp
colorbar
axis equal

if iT > 1
    caxis([cmin,cmax])
else
    cmin = min(min(W))-0.1*min(min(W));
    cmax = max(max(W))+0.1*max(max(W));
    
    caxis([cmin,cmax])
end

end % post : plotContour

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,W,cmin,cmax] = plotContourBox(o, variable, iT, xmin, xmax, ymin, ymax, ...
                        Nx, Ny, replace_inside, cmin, cmax, exact_sol)
    
[X, Y] = meshgrid(linspace(xmin, xmax, Nx), linspace(ymin, ymax, Ny));
[X,Y,W,cmin,cmax] = o.plotContour(variable, iT, X,Y, replace_inside, ...
                        cmin, cmax, exact_sol);

end % post : plotContourBox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,W,cmin,cmax] = plotContourCircle(o, variable, iT, rmin, rmax, cx, cy, ...
                        Nr, Nt, replace_inside, cmin, cmax, exact_sol)
   
theta = linspace(0, 2*pi, Nt);
r = linspace(rmin, rmax, Nr);
[R, T] = meshgrid(r, theta);

X = R.*cos(T) + cx;
Y = R.*sin(T) + cy;

[X,Y,W,cmin,cmax] = o.plotContour(variable, iT, X, Y, replace_inside,...
    cmin, cmax, exact_sol);

end % post : plotContourCircle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sxx, sxy, syx, syy] = computeStressOnWalls(o,iT)
    
oc = curve;
xWalls = oc.createWalls(o.prams.Nw, o.options);
walls = capsules([],xWalls);

RS = zeros(3*(o.prams.nw),1);
for k = 1:o.prams.nw-1
    RS(3*(k-1)+1:3*(k-1)+2) = o.stokeslets(2*k-1:2*k);
    RS(3*(k-1)+3) = o.rotlets(k);
end
% compute stresses from walls and Rotlets/Stokeslets
[stress1w, stress2w] = walls.stressTensor(o.eta_w(:,:,iT), RS, walls, true);

[sxxw,sxyw] = oc.getXY(stress1w);
[syxw,syyw] = oc.getXY(stress2w);

if o.prams.np > 0
    % compute stresses from fibers
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
    [stress1f, stress2f] = geom.stressTensor(o.eta_p(:,:,iT), [], walls, false);

    [sxxf,sxyf] = oc.getXY(stress1f);
    [syxf,syyf] = oc.getXY(stress2f);

else
    sxxf = zeros(o.prams.Nw,o.prams.nw);
    sxyf = zeros(o.prams.Nw,o.prams.nw);
    syxf = zeros(o.prams.Nw,o.prams.nw);
    syyf = zeros(o.prams.Nw,o.prams.nw);
end

sxx = sxxf+sxxw;
sxy = sxyf+sxyw;
syx = syxf+syxw;
syy = syyf+syyw;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [forceN, forceT] = computeWallForceComponents(o, iT)
    
forceN = zeros(1,o.prams.nw);
forceT = zeros(1,o.prams.nw);

[sxx, sxy, syx, syy] = computeStressOnWalls(o,iT);

oc = curve;

xWalls = oc.createWalls(o.prams.Nw, o.options);
walls = capsules([],xWalls);

[tx,ty] = oc.getXY(walls.xt);
[nx,ny] = oc.getXYperp(walls.xt);
sa = walls.sa;

for i = 1:o.prams.nw
    
    fdotn = sum(((sxx(:,i).*nx(:,i) + sxy(:,i).*ny(:,i)).*nx(:,i) + ...
                        (syx(:,i).*nx(:,i) + syy(:,i).*ny(:,i)).*ny(:,i)).*sa(:,i));
    
    fdott = sum(((sxx(:,i).*nx(:,i) + sxy(:,i).*ny(:,i)).*tx(:,i) + ...
                        (syx(:,i).*nx(:,i) + syy(:,i).*ny(:,i)).*ty(:,i)).*sa(:,i));
                    
    forceN(i) = fdotn*2*pi/o.prams.Nw;    
    forceT(i) = fdott*2*pi/o.prams.Nw;
end
    

end % post : computeWallForceComponents

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [force, torque] = computeWallTorque(o,iT)

force = zeros(2,o.prams.nw);
torque = zeros(1,o.prams.nw);

[sxx, sxy, syx, syy] = computeStressOnWalls(o,iT);
oc = curve;

xWalls = oc.createWalls(o.prams.Nw, o.options);
walls = capsules([],xWalls);

[x,y] = oc.getXY(walls.X);
[nx,ny] = oc.getXYperp(walls.xt);
sa = walls.sa;

for i = 1:o.prams.nw
    
    force(1,i) = sum((sxx(:,i).*nx(:,i) + sxy(:,i).*ny(:,i)).*sa(:,i));
    force(2,i) = sum((syx(:,i).*nx(:,i) + syy(:,i).*ny(:,i)).*sa(:,i));
    
    torque(i) = sum((y(:,i).*(sxx(:,i).*nx(:,i) + sxy(:,i).*ny(:,i)) - ...
                        x(:,i).*(syx(:,i).*nx(:,i) + syy(:,i).*ny(:,i))).*sa(:,i));
end
    
end % post : computeWallTorque
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stress_xx, stress_xy, stress_yx, stress_yy] = computeWallStress(o, iT)

oc = curve;
xWalls = oc.createWalls(o.prams.Nw, o.options);
walls = capsules([],xWalls);

sa = walls.sa;

[sxx, sxy, syx, syy] = computeStressOnWalls(o,iT);

for i = 2
    
    stress_xx = sum(sxx(:,i).*sa(:,i));
    stress_xy = sum(sxy(:,i).*sa(:,i));
    stress_yx = sum(syx(:,i).*sa(:,i));
    stress_yy = sum(syy(:,i).*sa(:,i));
end
    

end % post : computeWallStress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drag = computeWallDrag(o, iT)
    
drag = zeros(1,o.prams.nw);

oc = curve;
xWalls = oc.createWalls(o.prams.Nw, o.options);
walls = capsules([],xWalls);

[nx,ny] = oc.getXYperp(walls.xt);
sa = walls.sa;

[sxx, sxy,~,~] = computeStressOnWalls(o,iT);

for i = 1:o.prams.nw
    
    drag = sum((sxx(:,i).*nx(:,i) + sxy(:,i).*ny(:,i)).*sa(:,i));
    drag(i) = drag *2*pi/o.prams.Nw;    
end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stressTotal = evaluateVolumeAverageStress(o, iT)

stressTotal = zeros(2,2);
oc = curve;

if o.prams.np > 0
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));

    [stress1ff,stress2ff] = geom.stressTensor(o.eta_p(:,:,iT),[],geom,true);

    [x,y] = oc.getXY(geom.X);
    [nx,ny] = oc.getXYperp(geom.xt);
    sa = geom.sa;
    [sx1,sy1] = oc.getXY(stress1ff);
    [sx2,sy2] = oc.getXY(stress2ff);

    % compute stress from fibers on fibers
    for k = 1:o.prams.np
        ut = o.U(:,k,iT);
        omega = o.omega(k);
        xc = o.xc(:,k,iT);
        
        for i = 1:o.prams.Np
            
            rperp = [y(i,k)-xc(2);-(x(i,k)-xc(1))];
            u = ut + omega*rperp;
            
            f = [sx1(i,k), sx2(i,k); sy1(i,k), sy2(i,k)]*[nx(i,k);ny(i,k)];
            xof = [x(i,k);y(i,k)]*f';
            fox = f*[x(i,k),y(i,k)];
            
            uon = u*[nx(i,k),ny(i,k)];
            nou = [nx(i,k);ny(i,k)]*u';
            
            stressTotal = stressTotal + 0.5*(xof + fox - uon - nou)*sa(i,k)*2*pi/o.prams.Np;
        end
    end
end


if o.options.confined
    xWalls = oc.createWalls(o.prams.Nw, o.options);
    walls = capsules([],xWalls);
    
    RS = zeros(3*(o.prams.nw),1);
    for k = 1:o.prams.nw-1
        RS(3*(k-1)+1:3*(k-1)+2) = o.stokeslets(2*k-1:2*k);
        RS(3*(k-1)+3) = o.rotlets(k);
    end
    
    if o.prams.np > 0
        % stress from walls on fibers
        [stress1wf,stress2wf] = walls.stressTensor(o.eta_w(:,:,iT),RS,geom,false);
        [sx1,sy1] = oc.getXY(stress1wf);
        [sx2,sy2] = oc.getXY(stress2wf);
        
        for k = 1:o.prams.np
            for i = 1:o.prams.Np
                f = [sx1(i,k), sx2(i,k); sy1(i,k), sy2(i,k)]*[nx(i,k);ny(i,k)];
                xof = [x(i,k);y(i,k)]*f';
                
                stressTotal = stressTotal + xof*sa(i,k)*2*pi/o.prams.Np;
            end
        end
        
        % stress from fibers on walls
        [stress1fw,stress2fw] = geom.stressTensor(o.eta_p(:,:,iT),[],walls,false);
        
        [x,y] = oc.getXY(walls.X);
        [nx,ny] = oc.getXYperp(walls.xt);
        sa = walls.sa;  
        [sx1,sy1] = oc.getXY(stress1fw);
        [sx2,sy2] = oc.getXY(stress2fw);

        for k = 1:o.prams.nw
            for i = 1:o.prams.Nw

                f = [sx1(i,k), sx2(i,k); sy1(i,k), sy2(i,k)]*[nx(i,k);ny(i,k)];
                xof = [x(i,k);y(i,k)]*f';

                stressTotal = stressTotal + xof*sa(i,k)*2*pi/o.prams.NW;
            end
        end
    end
    
    [x,y] = oc.getXY(walls.X);
    [nx,ny] = oc.getXYperp(walls.xt);
    sa = walls.sa;  

    % stress from walls on walls
    [stress1ww,stress2ww] = walls.stressTensor(o.eta_w(:,:,iT),RS,walls,true);
    [sx1,sy1] = oc.getXY(stress1ww);
    [sx2,sy2] = oc.getXY(stress2ww);
    
    for k = 1:o.prams.nw        
        for i = 1:o.prams.NW
            
            f = [sx1(i,k), sx2(i,k); sy1(i,k), sy2(i,k)]*[nx(i,k);ny(i,k)];
            xof = [x(i,k);y(i,k)]*f';            
            
            stressTotal = stressTotal + xof*sa(i,k)*2*pi/o.prams.NW;
        end
    end
end

end % post : evaluateVolumeAverage

end %methods

end %classdef

