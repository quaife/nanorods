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
etaF
etaW
rotlets
stokeslets
wall_c

OUTPUTPATH_GIFS;

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = post(dataFile)
 
load(dataFile)

o.prams = prams;
o.options = options;
o.times = t;
o.xc = xc;
o.tau = tau;
o.U = U;
o.omega = omega;
o.etaF = etaF;
o.etaW = etaW;
o.rotlets = rot;
o.stokeslets = stokes;
%o.wall_c = wall_c;

o.OUTPUTPATH_GIFS = '../output/gifs/';

end %post : constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_fibres(o, iT, xmin, xmax, ymin, ymax)
    
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
    
    X = geom.getXY();
    oc = curve;
    [x, y] = oc.getXY(X);
    
    fill([x;x(1,:)],[y;y(1,:)],'k');

    axis equal

    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
    
    title(sprintf('t = %6.3f', o.times(iT)));
    
end % post : plot_fibres

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_fluid(o, iT, xmin, xmax, ymin, ymax, M, stride, bg_flow)
% plots the Stokes double layer potential over a meshgrid X, Y with M 
% points in each direction. 


[X, Y] = meshgrid(linspace(xmin, xmax, M), linspace(ymin, ymax, M));

geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
Utmp =  o.evaluateDLP(geom, o.etaF(:,:,iT), [X(:), Y(:)]);

% indentify points inside fiber
inside = geom.sortPts([X(:);Y(:)],o.options.ifmm);

%add in background flow
bg = bg_flow(X(:),Y(:));

Utmp(1:end/2) = Utmp(1:end/2) + bg(:,1); 
Utmp(end/2+1:end) = Utmp(end/2+1:end) + bg(:,2);
Utmp(inside==1) = 0;
U = reshape(Utmp(1:end/2), M, M);
V = reshape(Utmp(end/2+1:end), M, M);

quiver(X, Y, U, V, stride);

end % post : plot_fluid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,p] = plot_pressure(o,iT, xmin, xmax, ymin, ymax, M)


[X, Y] = meshgrid(linspace(xmin, xmax, M), linspace(ymin, ymax, M));
    
geomTar = capsules([],[X(:);Y(:)]);

oc = curve;
geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
p = geom.pressure(o.etaF(:,:,iT), [], geomTar);

% indentify points inside fiber
insideFibers = geom.sortPts([X(:);Y(:)],true);

if o.options.confined
    xWalls = oc.createWalls(o.prams.Nbd, o.options);

    walls = capsules(o.prams,xWalls);
    p = p + walls.pressure(o.etaW(:,:,iT), o.stokeslets(:,:,iT), geomTar);

    wallOuter = capsules(o.prams, xWalls(:,1));
    insideOuterWall = wallOuter.sortPts([X(:);Y(:)],true);
    
    % identify points outside outer wall
    p(insideOuterWall==0) = NaN;    
    
    if o.prams.Nbd > 1
        wallsInner = capsules(o.prams, xWalls(:,2:end));
        
        insideInnerWalls = wallsInner.sortPts([X(:);Y(:)],true);  
        % identify points inside inner walls
        p(insideInnerWalls==1) = NaN;
        X(insideInnerWalls==1) = NaN;
        Y(insideInnerWalls==1) = NaN;
    end
    
    X(insideOuterWall==0) = NaN;
    Y(insideOuterWall==0) = NaN;
end

p(insideFibers==1) = NaN;
X(insideFibers==1) = NaN;
Y(insideFibers==1) = NaN;

p = reshape(p,M,M);

surf(X,Y,p);
view(2);
shading interp
hold on

colorbar
caxis([-5,5]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_walls(o,iT)
    
    oc = curve;
    xWalls = oc.createWalls(o.prams.Nbd, o.options);
    
    hold on;
    for i = 1:o.prams.nbd
        plot([xWalls(1:end/2,i);xWalls(1,i)], ...
            [xWalls(end/2+1:end,i);xWalls(end/2+1,i)],'r', 'linewidth', 2);
    end

    ptsTrack = o.prams.tracker_fnc(o.times(iT));
    
    plot(ptsTrack(:,1),ptsTrack(:,2), 'b.', 'MarkerSize', 20); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = animated_gif(o, gif_options)
    
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
    xmin = min(min(o.xc(1,:,1:itmax))) - max(o.prams.lengths);
else
    if ~strcmp(gif_options.xmin, 'auto:frame')
        xmin = gif_options.xmin;
    end
end

if strcmp(gif_options.xmax, 'auto:all')
    xmax = max(max(o.xc(1,:,1:itmax))) + max(o.prams.lengths);
else
    if ~strcmp(gif_options.xmax, 'auto:frame')
        xmax = gif_options.xmax;
    end
end

if strcmp(gif_options.ymin, 'auto:all')
    ymin = min(min(o.xc(2,:,1:itmax))) - max(o.prams.lengths);
else
    if ~strcmp(gif_options.ymin, 'auto:frame')
        ymin = gif_options.ymin;
    end
end

if strcmp(gif_options.ymax, 'auto:all')
    ymax = max(max(o.xc(2,:,1:itmax))) + max(o.prams.lengths);
else
    if ~strcmp(gif_options.ymax, 'auto:frame')
        ymax = gif_options.ymax;
    end
end
   
% loop over all time steps
for i = 1:gif_options.stride:itmax
    
    clf;
    
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
      
    if gif_options.plot_pressure
        o.plot_pressure(i, xmin, xmax, ymin, ymax, gif_options.grid_pts);
    end
    
    if gif_options.plot_fluid
         o.plot_fluid(i, xmin, xmax, ymin, ymax, gif_options.grid_pts, ...
                        gif_options.velocity_grid_stride, gif_options.bg_flow); 
    end

    if o.options.confined
        o.plot_walls(i);
    end
    
    o.plot_fibres(i, xmin, xmax, ymin, ymax);
    
    drawnow;
    
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
        
        matlab2tikz([o.OUTPUTPATH_GIFS, gif_options.file_name, '/', ...
                 gif_options.file_name, '-', sprintf('%03d', i), '.tikz'],...
                 'height', '10cm', 'width', '12cm', 'standalone', true);
        
    else if strcmp(gif_options.file_type, 'gif')
            if i == 1;
                
                imwrite(imind, cm, ...
                    [o.OUTPUTPATH_GIFS,  gif_options.file_name, '.gif'], ...
                    'gif', 'Loopcount',inf, 'DelayTime', 0);
                
            else
                
                imwrite(imind, cm,...
                    [o.OUTPUTPATH_GIFS,  gif_options.file_name, '.gif'], ...
                    'gif','WriteMode','append', 'DelayTime',0);               
                
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
    NearStruct, @pot.exactStokesDL, @pot.exactStokesDL, geomTar,false,false);

end % post : evaluateDLP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stressTotal = evaluateTotalStress(o, iT)
   
stressTotal = zeros(2,2);
geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));

% [stress1ff,stress2ff] = geom.stressTensor(o.etaF,[],geom,true);
% 
oc = curve;
% [x,y] = oc.getXY(geom.X);
% [nx,ny] = oc.getXYperp(geom.xt);
% sa = geom.sa;
% [sx1,sy1] = oc.getXY(stress1ff);
% [sx2,sy2] = oc.getXY(stress2ff);

% compute stress from fibers on fibers
% for k = 1:o.prams.nv
% %    omega = o.omega(iT,k);
% %    u_trans = o.U(:,k,iT);
% %    xc = o.xc(1,k,iT);
% %    yc = o.xc(2,k,iT);
% %    
% %    u = u_trans + omega*[y(k) - yc; xc - x(k)];
% %    
%    for i = 1:o.prams.N
%        nox = [nx(i,k);ny(i,k)]*[x(i,k),y(i,k)];
% %        uon = u*[nx(i,k),ny(i,k)];
% %        nou = [nx(i,k);ny(i,k)]*u';
%        
%        stressTotal = stressTotal +...
%                 ([sx1(i,k), sx2(i,k); sy1(i,k), sy2(i,k)]*nox)*sa(i,k)*2*pi/o.prams.N;
%    end
% end



if o.options.confined
    xWalls = oc.createWalls(o.prams.Nbd, o.options);
    walls = capsules([],xWalls);
    
    RS = zeros(3*(o.prams.nbd),1);
    for k = 1:o.prams.nbd-1
        RS(3*(k-1)+1:3*(k-1)+2) = o.stokeslets(2*k-1:2*k);
        RS(3*(k-1)+3) = o.rotlets(k);
    end
    
%     % stress from walls on fibers  
%     [stress1wf,stress2wf] = walls.stressTensor(o.etaW,RS,geom,false);
%     [sx1,sy1] = oc.getXY(stress1wf);
%     [sx2,sy2] = oc.getXY(stress2wf);
%     for k = 1:o.prams.nv
%         for i = 1:o.prams.N
%             nox = [nx(i,k);ny(i,k)]*[x(i,k),y(i,k)];
%             
%             stressTotal = stressTotal +...
%                 2*pi/o.prams.N*[sx1(i,k), sx2(i,k); sy1(i,k), sy2(i,k)]*nox*sa(i,k);
%         end
%     end
    
    % stress from fibers on walls
    [x,y] = oc.getXY(walls.X);
    [nx,ny] = oc.getXYperp(walls.xt);
    sa = walls.sa;  
    
    [stress1fw,stress2fw] = geom.stressTensor(o.etaF,[],walls,false);
    [sx1,sy1] = oc.getXY(stress1fw);
    [sx2,sy2] = oc.getXY(stress2fw);
    
    for k = 1:o.prams.nbd
        for i = 1:o.prams.Nbd
  
            f = [sx1(i,k), sx2(i,k); sy1(i,k), sy2(i,k)]*[nx(i,k);ny(i,k)];
            xof = [x(i,k);y(i,k)]*f';
            
            stressTotal = stressTotal + xof*sa(i,k)*2*pi/o.prams.Nbd;
        end
    end
    
    % stress from walls on walls
    [stress1ww,stress2ww] = walls.stressTensor(o.etaW,RS,walls,true);
    [sx1,sy1] = oc.getXY(stress1ww);
    [sx2,sy2] = oc.getXY(stress2ww);
    
    for k = 1:o.prams.nbd
        for i = 1:o.prams.Nbd
            nox = [nx(i,k);ny(i,k)]*[x(i,k),y(i,k)];
            
            stressTotal = stressTotal +...
                2*pi/o.prams.Nbd*[sx1(i,k), sx2(i,k); sy1(i,k), sy2(i,k)]*nox*sa(i,k);
        end
    end
end


end % post : evaluateTotalStress

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[stress1, stress2] = evaluateStress(o, iT, X)
% evaluates the stress tensor at time step iT at target points X, given a
% boundary type, either 'fibers', or 'walls'

geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
geomTar = capsules([],X);

RS = [];
[stress1,stress2] = geom.stressTensor(o.etaF,RS,geomTar,false);

if o.options.confined

    oc = curve;
    xWalls = oc.createWalls(o.prams.Nbd, o.options);
    walls = capsules([],xWalls);
    [stress1w,stress2w] = walls.stressTensor(o.etaW,RS,geomTar,false);
    
    stress1 = stress1 + stress1w;
    stress2 = stress2 + stress2w;
end

end % post : evaluateDLP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pressure = evaluatePressure(o, iT, X, boundary_type)
% evaluates the stress tensor at time step iT at target points X, given a 
% boundary type, either 'fibers', or 'walls'

switch boundary_type
    case 'fibers'
        geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
        eta = o.etaF;
        RS = [];
        
    case 'walls'
        
end

geomTar = capsules([],X);

pressure = geom.pressure(eta,RS,geomTar);

end % post : evaluateDLP
end %methods

end %classdef

