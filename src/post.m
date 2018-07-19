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
force_p
torque_p

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

if exist('force_p')
    o.force_p = force_p;
    o.torque_p = torque_p;
end
%o.wall_c = wall_c;

o.OUTPUTPATH_GIFS = '../output/gifs/';

end %post : constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotFibres(o, iT)
    
        geom0 = capsules(o.prams, o.xc(:,:,1), o.tau(1,:));
        fill(geom0.X(1:end/2,:), geom0.X(end/2+1:end,:), 'b');
    
%     ind = [20,10];
%     theta = -(0:o.prams.Np-1)'*2*pi/o.prams.Np;
%     X = zeros(2*o.prams.Np, o.prams.np);
%     for k = 1:o.prams.np
%         x = 0.5*o.prams.lengths*(1+0.02*cos(ind(k)*theta)).*cos(theta); y = 0.5*o.prams.widths*(1+0.02*cos(ind(k)*theta)).*sin(theta);
%         X(:,k) = [x*cos(o.tau(iT,k)) - y*sin(o.tau(iT,k)) + o.xc(1,k,iT);
%             x*sin(o.tau(iT,k)) + y*cos(o.tau(iT,k)) + o.xc(2,k,iT)];
%     end
%     
%     geom = capsules(o.prams, X);
    
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
    fill3([geom.X(1:end/2,:);geom.X(1,:)], ...
            [geom.X(end/2+1:end,:);geom.X(end/2+1,:)],...
            1*10000000*ones(o.prams.Np+1, 1),'k');
    %plot([geom.X(1:end/2,:); geom.X(1,:)],[geom.X(end/2+1:end,:);geom.X(end/2+1,:)], 'b');
    
%     
	 x = -12*ones(4,1);
     y = [-2, -1, 1, 2];
	 v = zeros(4,1);
     u = [-2; -1; 1; 2];
    
     if o.times(iT) > 10
         u = -u;
     end
     
     quiver(x,y,u,v)
     
%     hold on;
%     counter_particles = find(o.torque_p(:,1) < 0);
%     for i = 1:geom.n
%         if max(counter_particles == i)
%             color = 'b';
%         else
%             color = 'k';
%         end
%         
%         h = fill3([geom.X(1:end/2,i);geom.X(1,i)], ...
%             [geom.X(end/2+1:end,i);geom.X(end/2+1,i)],...
%             1*10000000*ones(o.prams.Np+1, 1),color);
%         
%         %set(h, 'facealpha', 0.25);
%         %plot(geom.X(1,i), geom.X(end/2+1,i), 'mo');
%     end

    contact_particles = find(max(abs(o.force_p(:,:,iT))) > 1e-8);
    if ~isempty(contact_particles)
        fill3([geom.X(1:end/2,contact_particles);geom.X(1,contact_particles)], ...
            [geom.X(end/2+1:end,contact_particles);geom.X(end/2+1,contact_particles)],...
            1*10000000*ones(o.prams.Np+1, 1),'g');
    end
%     plot([geom.X(1:end/2,contact_particles); geom.X(1,contact_particles)],...
%            [geom.X(end/2+1:end,contact_particles); geom.X(end/2+1,contact_particles)], 'g', 'linewidth', 4);
    
    
end % post : plot_fibres

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotParticleTracers(o, iT)
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
    hold on
    T = 1:iT;
    %for i = 1%:geom.n
       plot3(reshape(o.xc(1,1,T),1,length(T)), reshape(o.xc(2,1,T),1,length(T)), 20000000*ones(length(T)), 'b', 'linewidth', 2); 
       plot3(reshape(o.xc(1,2,T),1,length(T)), reshape(o.xc(2,2,T),1,length(T)), 20000000*ones(length(T)), 'c', 'linewidth', 2); 
      % plot(reshape(o.xc(1,3,T),1,length(T)), reshape(o.xc(2,3,T),1,length(T)), 'm', 'linewidth', 2); 
%        plot(reshape(o.xc(1,4,T),1,length(T)), reshape(o.xc(2,4,T),1,length(T)), 'g', 'linewidth', 2); 
    %end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotWalls(o,iT)
    
    oc = curve;
    xWalls = oc.createWalls(o.prams.Nw, o.options);
    walls = capsules(o.prams, xWalls);
    
    hold on;
    for i = 1:o.prams.nw
        plot3([xWalls(1:end/2,i);xWalls(1,i)], ...
            [xWalls(end/2+1:end,i);xWalls(end/2+1,i)],...
            100*ones(o.prams.Nw+1, 1),'r', 'linewidth', 4);
    end
    

    for i = 3:walls.n
        fill3([walls.X(1:end/2,i);walls.X(1,i)], ...
            [walls.X(end/2+1:end,i);walls.X(end/2+1,i)],...
            10000000*ones(o.prams.Nw+1, 1),'k')
    end
    
    if ~isempty(o.prams.tracker_fnc)
        ptsTrack = o.prams.tracker_fnc(o.times(iT));
        
        plot3(ptsTrack(:,1),ptsTrack(:,2), 1000*ones(size(ptsTrack,1),1),...
            'b.', 'MarkerSize', 20);
    end

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
   
if ~isempty(gif_options.dt)
    irange = 1;
    t = o.times(1);
    while t < o.times(end)
        t = t + gif_options.dt;
        t_target = t;
        [~,irange(end+1)] = min(abs(o.times - t_target));
    end
else
    irange = 1:gif_options.stride:itmax;
end

% loop over all time steps
for i = irange
    
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
            o.plotContourCircle(gif_options.contour_field, i, 0.5, 7, 0, 0, 8, 5, ...
                                true, cmin, cmax, []);
        else
            [~,~,~,cmin,cmax] = o.plotContourCircle(gif_options.contour_field, i, 0.5, 7, ...
                                    0, 0, 8, 5, true, [], [], []);
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
    
    if gif_options.tracers
        o.plotParticleTracers(i);
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
                 gif_options.file_name, '-', sprintf('%04d', i), '.tikz'],...
                 'height', '10cm', 'width', '12cm', 'standalone', true, 'floatFormat', '%.3f');
        
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
function animatedGifExtended(o, gif_options)
    
h = figure();
set(h,'Units','normalized');
set(h,'Position',[0 0 1 1]);

[~,~,m] = size(o.xc);

% get axis limits
if strcmp(gif_options.xmin, 'auto:all')
    xmin = min(min(o.xc(1,:,:))) - o.prams.lengths;
else
    if ~strcmp(gif_options.xmin, 'auto:frame')
        xmin = gif_options.xmin;
    end
end

if strcmp(gif_options.xmax, 'auto:all')
    xmax = max(max(o.xc(1,:,:))) + o.prams.lengths;
else
    if ~strcmp(gif_options.xmax, 'auto:frame')
        xmax = gif_options.xmax;
    end
end

if strcmp(gif_options.ymin, 'auto:all')
    ymin = min(min(o.xc(2,:,:))) - o.prams.lengths;
else
    if ~strcmp(gif_options.ymin, 'auto:frame')
        ymin = gif_options.ymin;
    end
end

if strcmp(gif_options.ymax, 'auto:all')
    ymax = max(max(o.xc(2,:,:))) + o.prams.lengths;
else
    if ~strcmp(gif_options.ymax, 'auto:frame')
        ymax = gif_options.ymax;
    end
end
if strcmp(gif_options.itmax, 'all')
    itmax = m;
else
    itmax = gif_options.itmax;
end

% loop over all time steps
for i = 1:gif_options.stride:itmax-1
    
        
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
    

    subplot(3,2,1);
    o.plotFibres(i);
    axis equal
    title(sprintf('t = %6.3f', o.times(i)));
    
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    
    if ~gif_options.axis 
        axis off
    end
    
    view(2);
    
    subplot(3,2,2);
    plot(o.eta_p(:,:,i));
            
    subplot(3,2,3);
    z = fftshift(fft(o.eta_p(1:end/2,1,i) + 1i*o.eta_p(end/2+1:end,1,i)));
    semilogy(abs(z), 'b');
    ylim([0,1e5])
    
    subplot(3,2,4);
    z = fftshift(fft(o.eta_p(1:end/2,2,i) + 1i*o.eta_p(end/2+1:end,2,i)));
    semilogy(abs(z), 'r');
    ylim([0,1e5])
    
    subplot(3,2,5);
    plot(o.force_p(:,:,i))

    subplot(3,2,6);
    plot(o.torque_p(:,i));

    
    drawnow
    frame = getframe(h);
    im = frame2im(frame);
    
    [imind, cm] = rgb2ind(im,256);
    
    if i == 1;
                
        imwrite(imind, cm, ...
            [o.OUTPUTPATH_GIFS,  gif_options.file_name, '.gif'], ...
            'gif', 'Loopcount',inf, 'DelayTime', 0.3);

    else

        imwrite(imind, cm,...
            [o.OUTPUTPATH_GIFS,  gif_options.file_name, '.gif'], ...
            'gif','WriteMode','append', 'DelayTime',0.3);               

    end
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = evaluateDLP(o, geom, eta, X)
% evaluates the Stokes double layer potential at target points X, given a 
% geometry geom and a density function eta


geomTar = capsules([],X);

pot = poten(geom.N);
D = pot.stokesDLmatrix(geom);

if o.options.near_singular
    
    [~,NearStruct] = geom.getZone(geomTar,2);

    D = pot.stokesDLmatrix(geom);
    DLP = @(X) pot.exactStokesDLdiag(geom,D,X) - 1/2*X;
    u = pot.nearSingInt(geom, eta, DLP,[],...
        NearStruct, @pot.exactStokesDLfmm, @pot.exactStokesDL, geomTar,false,false);
else
    [~, u] = pot.exactStokesDL(geom, eta, D, X, 1);
end

end % post : evaluateDLP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = evaluateRotStokes(o, iT, X)
    
% evalute rotlets and Stokeslets from solid walls
xi = o.rotlets(iT,:);
lambda = o.stokeslets(:,:,iT);

oc = curve;
[xWalls,~] = oc.createWalls(o.prams.Nw, o.options);
walls = capsules([],xWalls);

[x,y] = oc.getXY(X);

u = zeros(2,length(x));

for i = 1:length(xi)
    r = [x - walls.center(1,i+1); y - walls.center(2,i+1)];

    for j = 1:length(x)

        rho = norm(r(:,j));
        u(:,j) = u(:,j) + xi(i)*[r(2,j);-r(1,j)]/rho^2/(4*pi);

        ror = r(:,j)*r(:,j)';

        u(:,j) = u(:,j) + (ror/rho^2 - log(rho)*eye(2))*lambda(:,i)/(4*pi);
    end
end

% evalute rotlets and Stokeslets from particles
xi = o.torque_p(:,iT);
lambda = o.force_p(:,:,iT);
cen = o.xc(:,:,iT);

for i = 1:length(xi)
    r = [x - cen(1,i); y - cen(2,i)];
    for j = 1:length(x)

        rho = norm(r(:,j));
        u(:,j) = u(:,j) + xi(i)*[r(2,j);-r(1,j)]/rho^2/(4*pi); 
        ror = r(:,j)*r(:,j)';

        u(:,j) = u(:,j) + (ror/rho^2 - log(rho)*eye(2))*lambda(:,i)/(4*pi);
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
function [stress1, stress2] = evaluateStress(o, iT, X)
% evaluates the stress tensor at time step iT at target points X, given a
% boundary type, either 'fibers', or 'walls'

if o.prams.np > 0
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
else
    geom = capsules(o.prams, [], []);
end

geomTar = capsules([],X);

RS = zeros(3*(o.prams.np),1);
for k = 1:o.prams.np
    RS(3*k+1:3*k+2) = o.force_p(:,k,iT);
    RS(3*k+3) = o.torque_p(k,iT);
end

if geom.n > 0
    [stress1,stress2] = geom.stressTensor(o.eta_p(:,:,iT),RS,geomTar,false);
else
    stress1 = zeros(2,geomTar.n);
    stress2 = zeros(2,geomTar.n);
end

if o.options.confined

    RS = zeros(3*(o.prams.nw),1);
    for k = 2:o.prams.nw
        RS(3*k+1:3*k+2) = o.stokeslets(:,k-1,iT);
        RS(3*k+3) = o.rotlets(iT,k-1);
    end
    
    oc = curve;
    [xWalls,~] = oc.createWalls(o.prams.Nw, o.options);
    walls = capsules([],xWalls);
    
    RS = [];
    
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

end % post : velocityQuiverBox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,U,V] = velocityQuiverPolar(o, iT, rmin, rmax, ...
                        Nr, Ntheta, stride, bg_flow)
    
[R, THETA] = meshgrid(linspace(rmin, rmax, Nr), linspace(0, 2*pi, Ntheta));

X = R.*cos(THETA);
Y = R.*sin(THETA);

[X,Y,U,V] = o.velocityQuiver(iT, X, Y, stride, bg_flow);

o.plotWalls(iT)
axis equal


end % post : velocityQuiverBox

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

Xnodes = linspace(min(min(X)), max(max(X)), 100);
Ynodes = linspace(min(min(Y)), max(max(Y)), 100);
[zgrid,xgrid,ygrid] = gridfit(X(:),Y(:),W(:),Xnodes,Ynodes);

% remove data inside r < 0.5
for i = 1:100
    for j = 1:100
        if xgrid(i,j)^2 + ygrid(i,j)^2 < 0.4^2
            zgrid(i,j) = nan;
        end
    end
end

%W = reshape(W,M,N);

% if min(size(X))==1
%     nanflags = isnan(W);
%     X(nanflags) = [];
%     Y(nanflags) = [];
%     W(nanflags) = [];
%     
%     F = scatteredInterpolant(X',Y', W');
%     F.Method = 'natural';
%     X1 = linspace(min(X),max(X),100);
%     Y1 = linspace(min(Y),max(Y),100);
%     [X,Y] = meshgrid(X1,Y1);
%     
%     W = F(X,Y);
% end

%surf(X,Y,log10(abs(W)));
surf(xgrid, ygrid, log10(abs(zgrid)));
view(2);
pause(0.1);
shading interp
%colorbar
axis equal

cmin = 2;
cmax = 8;
% cmin = min(min(zgrid))-0.1*min(min(zgrid));
% cmax = max(max(zgrid))+0.1*max(max(zgrid));

if ~isempty(cmin)
    caxis([cmin,cmax])
else
    cmin = min(min(log10(abs(W))))-0.1*min(min(log10(abs(W))));
    cmax = max(max(log10(abs(W))))+0.1*max(max(log10(abs(W))));
    
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
r_possible = linspace(rmin, rmax, Nr);

theta = linspace(0,2*pi,Nt);
r(1:Nt) = r_possible(1);
spacing = 2*pi*r(1)^2/Nt;

for i = 2 : length(r_possible)
    N = ceil(2*pi*r_possible(i)/spacing);
    
    theta_new = linspace(0, 2*pi, N);
    for j = 1:length(theta_new)
        theta(end+1) = theta_new(j);
        r(end+1) = r_possible(i);
    end
end

R = r;
T = theta;


%[R, T] = meshgrid(r_possible, theta);

X = R.*cos(T) + cx;
Y = R.*sin(T) + cy;

[X,Y,W,cmin,cmax] = o.plotContour(variable, iT, X, Y, replace_inside,...
    cmin, cmax, exact_sol);

o.plotWalls(iT);
axis equal

hold on
o.plotFibres(iT)

xlim([-6,6]);
ylim([-6,6]);

axis off

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
function K = plot_kernel(o, iT, x_target)
    
    oc = curve;
    xWalls = oc.createWalls(o.prams.Nw, o.options);
    walls = capsules([],xWalls);
    [nx,ny] = oc.getXYperp(walls.xt);
    [x,y] = oc.getXY(walls.X);
    [etax, etay] = oc.getXY(o.eta_w(:,iT,1));
    
    rx = x - x_target(1);
    ry = y - x_target(2);
    r = [rx, ry];
    K = zeros(walls.N, 2);
    
    for i = 1:walls.N
       
        ror = r(i,:)'*r(i,:);
        rin = r(i,:)*[nx(i);ny(i)];
        
        K(i,:) = rin*ror*[etax(i); etay(i)]/(r(i,1)^2 + r(i,2)^2).^2;
    
    end

    clf;
    t = linspace(0,2*pi,walls.N);
    plot(t, K);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stressFF, stressWF] = evaluateVolumeAverageStress(o, iT)


stressFF = zeros(2,2);
stressFW = zeros(2,2);
stressWF = zeros(2,2);
stressWW = zeros(2,2);

oc = curve;

if o.prams.np > 0
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
    
    RS = zeros(3*(o.prams.np),1);
    for k = 1:o.prams.np
        RS(3*(k-1)+1:3*(k-1)+2) = o.force_p(:,k,iT);
        RS(3*(k-1)+3) = o.torque_p(k,iT);
    end
    
    %RS = [];
    [stress1ff,stress2ff] = geom.stressTensor(o.eta_p(:,:,iT),RS,geom,true);
    
    [x,y] = oc.getXY(geom.X);
    [nx,ny] = oc.getXYperp(geom.xt);
    sa = geom.sa;
    [sx1,sy1] = oc.getXY(stress1ff);
    [sx2,sy2] = oc.getXY(stress2ff);
    
    % compute stress from fibers on fibers
    for k = 1:o.prams.np
        
        for i = 1:o.prams.Np

            f = [sx1(i,k), sx2(i,k); sy1(i,k), sy2(i,k)]*[nx(i,k);ny(i,k)];
            xof = [x(i,k);y(i,k)]*f';
            
            t = atan2(y(i,k),x(i,k));
            Q = [cos(t), sin(t); -sin(t), cos(t)];
            %compute stress tensor in polar coordinates
            stressFF = stressFF + Q*xof*Q'*sa(i,k)*2*pi/o.prams.Np;
        end
    end
end

if o.options.confined
    xWalls = oc.createWalls(o.prams.Nw, o.options);
    walls = capsules([],xWalls);
    
    RS = zeros(3*(o.prams.nw),1);
    for k = 2:o.prams.nw
        RS(3*(k-1)+1:3*(k-1)+2) = o.stokeslets(:,k-1,iT);
        RS(3*(k-1)+3) = o.rotlets(iT,k-1);
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
                t = atan2(y(i,k),x(i,k));
                Q = [cos(t), sin(t); -sin(t), cos(t)];
            
                stressWF = stressWF + Q*xof*Q'*sa(i,k)*2*pi/o.prams.Np;
            end
        end
        
        %         % stress from fibers on walls
        %         [stress1fw,stress2fw] = geom.stressTensor(o.eta_p(:,:,iT),[],walls,false);
        %
        %         [x,y] = oc.getXY(walls.X);
        %         [nx,ny] = oc.getXYperp(walls.xt);
        %         sa = walls.sa;
        %         [sx1,sy1] = oc.getXY(stress1fw);
        %         [sx2,sy2] = oc.getXY(stress2fw);
        %
        %         for k = 1:o.prams.nw
        %             for i = 1:o.prams.Nw
        %
        %                 f = [sx1(i,k), sx2(i,k); sy1(i,k), sy2(i,k)]*[nx(i,k);ny(i,k)];
        %                 fox = f*[x(i,k),y(i,k)];
        %
        %                 stressFW = stressFW + fox*sa(i,k)*2*pi/o.prams.Nw;
        %             end
        %         end
    end
    %
%     [x,y] = oc.getXY(walls.X);
%     [nx,ny] = oc.getXYperp(walls.xt);
%     sa = walls.sa;  
% 
%     % stress from walls on walls
%     [stress1ww,stress2ww] = walls.stressTensor(o.eta_w(:,:,iT),RS,walls,true);
%     [sx1,sy1] = oc.getXY(stress1ww);
%     [sx2,sy2] = oc.getXY(stress2ww);
%     
%     for k = 1:o.prams.nw        
%         for i = 1:o.prams.Nw
%             
%             f = [sx1(i,k), sx2(i,k); sy1(i,k), sy2(i,k)]*[nx(i,k);ny(i,k)];
%             fox = f*[x(i,k),y(i,k)];            
%             
%             stressWW = stressWW + fox*sa(i,k)*2*pi/o.prams.Nw;
%         end
%     end
end

%stressTotal = stressFF + stressWF;

end % post : evaluateVolumeAverage

end %methods

end %classdef

