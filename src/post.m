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

EPS;
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

o.EPS = 0.1;
o.OUTPUTPATH_GIFS = '../output/gifs/';

end %post : constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_fibres_with_stats(o, iT, xmin, xmax, ymin, ymax)
    
    ax1 = axes('Position', [0 0 1 1], 'Visible', 'off');
    ax2 = axes('Position', [0.1 0.3 0.8 0.6]);

    
    geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
    X = geom.getXY();
    oc = curve;
    [x, y] = oc.getXY(X);
    
    axes(ax2)
    fill([x;x(1,:)],[y;y(1,:)],'k');
    
    axis equal
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
    
    title(sprintf('t = %6.3f', o.times(iT))); 
    axes(ax1)
    text(0.3, 0.15, {['N       : ' num2str(o.prams.N)], ['nv      : ' num2str(o.prams.nv)], ...
        ['order : ', num2str(o.prams.order(1))]});
    
    text(0.5, 0.15, {['timestep order : ' num2str(o.prams.tstep_order)], ...
                    ['FMM      : ' num2str(o.prams.fmm)], ...
                    ['near singular : ', num2str(o.prams.near_singular)]});       
                
end % post : plot_fibres_with_stats

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
function [] = plot_fluid(o, iT, xmin, xmax, ymin, ymax, epsilon)
% plots the Stokes double layer potential over a meshgrid X, Y. The DLP will
% only be evaluated if it is at least epsilon away from all fibers.
%
% TO DO : add different types of background flow, right now only
% extensional

%theta = (0:o.N-1)'*2*pi/o.N;

% for k = 1:o.nv
%    f(1:o.N,k) = 0*sin(theta);
%    f(o.N+1:2*o.N,k) = cos(theta);
% end

% f(1:o.N,1) = exp(sin(theta));
% f(o.N+1:2*o.N,1) = cos(theta);
% 
% % f(1:o.N,1) = ones(size(theta));
% % f(o.N+1:2*o.N,1) = ones(size(theta));
% 
% f(1:o.N,2) = exp(sin(theta));
% f(o.N+1:2*o.N,2) = exp(cos(theta));


% o.eta = f(:,:);

M = 30;
% X = linspace(xmin, xmax, M);
% Y = ymax*ones(size(X));
% 
% X = X(1:end-1);
% Y = Y(1:end-1);

[X, Y] = meshgrid(linspace(xmin, xmax, M), linspace(ymin, ymax, M));




geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
targetPts.nv = 1;
targetPts.N = M*M;
targetPts.X = [X(:), Y(:)];



Utmp =  o.evaluateDLP(geom, o.etaF(:,:,iT), X(:), Y(:));
%[~, nearOther] = geom.getZone(targetPts, 2);

% J = find(nearOther.zone{1}(:,2)==1);
% Utmp(J) = 0;
% Utmp(J+length(Utmp)/2) = 0;

U = reshape(Utmp(1:end/2), M, M);
V = reshape(Utmp(end/2+1:end), M, M);
U = U + 5*Y; %add in background flow

% U = Utmp(1:end/2);
% V = Utmp(end/2+1:end);

quiver(X, Y, U, V, 1);

%surf(X,Y,U);
% plot(X, U, '-o','linewidth', 2);

%view(2)

%shading interp
hold on 

o.plot_fibres(iT, xmin, xmax, ymin, ymax);
% 
% xlim([xmin, xmax]);
% ylim([ymin, ymax]);
%axis equal

title(sprintf('t = %6.3f', o.times(iT)));

% figure(2);
% hold on
% surf(X,Y,V);
% 
% % plot(X, V, '-o','linewidth', 2);
% view(2)
% 
% shading interp
% hold on
% 
% o.plot_fibres(iT, xmin, xmax, ymin, ymax);
% 
% % xlim([xmin, xmax]);
% % ylim([ymin, ymax]);
% %axis equal
% 
% title(sprintf('V at t = %6.3f', o.times(iT)));

end % post : plot_fluid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,p] = plot_pressure(o,iT)

M = 50;
N = 500;
r = linspace(5, 10, M);
theta = linspace(0,2*pi,N);

[R, T] = meshgrid(r, theta);

X = R.*cos(T);
Y = R.*sin(T);
geomTar = capsules([],[X(:);Y(:)]);

oc = curve;
geom = capsules(o.prams, o.xc(:,:,iT), o.tau(iT,:));
xWalls = oc.createWalls(o.prams.Nbd, o.options);

walls = capsules(o.prams,xWalls);

p = geom.pressure(o.etaF(:,:,iT), [], geomTar, true, 'DLP');
p = p + walls.pressure(o.etaW(:,:,iT), o.stokeslets(:,:,iT), geomTar, true, 'DLP');
p = reshape(p,N,M);

for i = 1:N
    for j = 1:M
    
        if any(sum((o.xc(:,:,iT) - repmat([X(i,j);Y(i,j)], 1, geom.nv)).^2) < (max(o.prams.widths/2))^2)
            p(i,j) = 0;
        end
    end
end

surf(X,Y,p);
view(2);
shading interp
hold on

[x, y] = oc.getXY(geom.X);
[n,m] = size([x;x(1,:)]);
fill3([x;x(1,:)],[y;y(1,:)],100*ones(n,m),'k');

colorbar
%caxis([-5,5]);
axis equal

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_walls(o,iT)
    
    oc = curve;
    xWalls = oc.createWalls(o.prams.Nbd, o.options);
    
    hold on;
    for i = 1:o.prams.nbd
        plot([xWalls(1:end/2,i);xWalls(1,i)], [xWalls(end/2+1:end,i);xWalls(end/2+1,i)],...
                'r', 'linewidth', 2);
    end

    ptsTrack = o.prams.tracker_fnc(o.times(iT));
    
    plot(ptsTrack(:,1),ptsTrack(:,2), 'b.', 'MarkerSize', 20); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = animated_gif(o, gif_options)
    
h = figure();
set(h,'Units','normalized');
set(h,'Position',[0 0 1 1]);

% get axis limits
if strcmp(gif_options.xmin, 'auto:all')
    xmin = min(min(o.xc(1,:,1:gif_options.itmax))) - max(o.prams.lengths);
else
    if ~strcmp(gif_options.xmin, 'auto:frame')
        xmin = gif_options.xmin;
    end
end

if strcmp(gif_options.xmax, 'auto:all')
    xmax = max(max(o.xc(1,:,1:gif_options.itmax))) + max(o.prams.lengths);
else
    if ~strcmp(gif_options.xmax, 'auto:frame')
        xmax = gif_options.xmax;
    end
end

if strcmp(gif_options.ymin, 'auto:all')
    ymin = min(min(o.xc(2,:,1:gif_options.itmax))) - max(o.prams.lengths);
else
    if ~strcmp(gif_options.ymin, 'auto:frame')
        ymin = gif_options.ymin;
    end
end

if strcmp(gif_options.ymax, 'auto:all')
    ymax = max(max(o.xc(2,:,1:gif_options.itmax))) + max(o.prams.lengths);
else
    if ~strcmp(gif_options.ymax, 'auto:frame')
        ymax = gif_options.ymax;
    end
end
   
% loop over all time steps
for i = 1:gif_options.stride:gif_options.itmax
    
    clf;
    
    if o.options.confined
        o.plot_walls(i);
    end
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
    
    if gif_options.plot_fluid
         o.plot_fluid(i, xmin, xmax, ymin, ymax, o.EPS); 
    else
        o.plot_fibres(i, xmin, xmax, ymin, ymax);
    end
    
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
        
        matlab2tikz([o.OUTPUTPATH_GIFS, gif_options.file_name, '/', gif_options.file_name, '-', sprintf('%03d', i), '.tikz'], ...
            'height', '10cm', 'width', '12cm', 'standalone', true);
        
    else if strcmp(gif_options.file_type, 'gif')
            if i == 1;
                
                imwrite(imind, cm, [o.OUTPUTPATH_GIFS,  gif_options.file_name, '.gif'], 'gif', ...
                    'Loopcount',inf, 'DelayTime', 0);
                
            else
                
                imwrite(imind, cm, [o.OUTPUTPATH_GIFS,  gif_options.file_name, '.gif'], 'gif',...
                    'WriteMode','append', 'DelayTime',0);               
                
            end
        end
    end
end

end % post : aimated_gif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stats = calculate_stats(o, fibers)
%provide statistics for fibers

stats.time = o.times;


stats.lengths = o.lengths(fibers);
stats.widths = o.widths(fibers);
stats.aspect_ratios = o.lengths(fibers)./o.widths(fibers);

stats.centre_x = o.centres_x(:,fibers);
stats.centre_y = o.centres_y(:,fibers);
stats.orientations = wrapTo2Pi(o.orientations(:,fibers));
stats.vel_trans_x = o.u_x(:,fibers);
stats.vel_trans_y = o.u_y(:,fibers);
stats.vel_ang = o.omega(:,fibers);
stats.speed = sqrt(stats.vel_trans_x.^2 + stats.vel_trans_y.^2);

%average quanities
stats.centre_x_avg = mean(stats.centre_x,2);
stats.centre_y_avg = mean(stats.centre_y,2);
stats.orientations_avg = mean(stats.orientations,2);
stats.vel_trans_x_avg = mean(stats.vel_trans_x,2);
stats.vel_trans_y_avg = mean(stats.vel_trans_y,2);
stats.vel_ang_avg = mean(stats.vel_ang,2);
stats.speed_avg = mean(stats.speed, 2);

%standard deviations
stats.centre_x_std = std(stats.centre_x, 0, 2);
stats.centre_y_std = std(stats.centre_y, 0, 2);
stats.orientations_std = std(stats.orientations, 0, 2);
stats.vel_trans_x_std = std(stats.vel_trans_x, 0, 2);
stats.vel_trans_y_std = std(stats.vel_trans_y, 0, 2);
stats.vel_ang_std = std(stats.vel_ang, 0, 2);
stats.speed_std = std(stats.speed, 0, 2);

%other quantities
stats.dist = sqrt((stats.centre_x(end,:)-stats.centre_x(1,:)).^2 -...
                (stats.centre_y(end,:) - stats.centre_y(1,:)).^2);
stats.n_rotations = (o.orientations(end,:) - o.orientations(1,:))/(2*pi);

%stats.orientations = 2*pi*rand(length(stats.time), 10000);
% for i = 1:length(stats.time)
% stats.orientations(i,1:10000) = linspace(0,2*pi,10000);
% end

%histogram and tensors
%apply symmetry around 0


stats.hist_thetas = linspace(0,2*pi,length(fibers));
%stats.hist_thetas = linspace(0,2*pi,1000);

stats.hist_p = hist(stats.orientations', stats.hist_thetas)';

for i = 1:floor(length(stats.hist_thetas)/2)
    ptmp = (stats.hist_p(:,i) + stats.hist_p(:,end-i+1))/2;
    stats.hist_p(:,i) = ptmp;
    stats.hist_p(:,end-i+1) = ptmp;
end

%normalize histogram to create "real" pdf
dtheta = stats.hist_thetas(2) - stats.hist_thetas(1);

stats.pdf_actual = stats.hist_p;
for i = 1:length(stats.time)
    stats.pdf_actual(i,:) = stats.hist_p(i,:)/(dtheta*sum(stats.hist_p(i,:)));
end

%calculate second order moment
stats.a2 = zeros(2,2,length(stats.time));
for i = 1:length(stats.time)
    for j = 1:length(stats.hist_thetas)
        ptmp = [cos(stats.hist_thetas(j)); sin(stats.hist_thetas(j))];
        stats.a2(:,:,i) = stats.a2(:,:,i) + dtheta*(ptmp*ptmp')*stats.pdf_actual(i,j);
    end
end

%calculate fourth order moment
stats.a4 = zeros(2,2,2,2,length(stats.time));

for i = 1:length(stats.time)
    for j = 1:length(stats.hist_thetas)
        ptmp = [cos(stats.hist_thetas(j)), sin(stats.hist_thetas(j))];
       
        %probably a better way to do this (at least more astetically
        %pleasing)
        for iT = 1:2
            for jT = 1:2
                for kT = 1:2
                    for nT = 1:2
                        stats.a4(iT,jT,kT,nT,i) = stats.a4(iT,jT,kT,nT,i) + ...
                            dtheta*ptmp(iT)*ptmp(jT)*ptmp(kT)*ptmp(nT);
                    end
                end
            end
        end
    end
end

end %post : stats

function prob = recover_probability(o, theta, a2, a4)
    %uses the second order tensor a2 calculated in the stats routine to
    %recover the probablity of finding a fiber with orientation angle
    %theta. See Advani and Tucker page 760
    
    p = [cos(theta);sin(theta)];
    p_outer = p*p';
    
    b2 = a2 - 0.5*eye(2);
    f2 = p_outer - 0.5*eye(2);
    prob = 1/(2*pi) + (2/pi)*trace(b2*f2');
    
    %if fourth order tensor given add on contribution
    if ~isempty(a4);        
        
        C1 = zeros(2,2,2,2);
        C1(1,1,:,:) = a2;
        C1(2,2,:,:) = a2;
        
        C2 = zeros(2,2,2,2);
        C2(1,:,1,:) = a2;
        C2(2,:,2,:) = a2;
        
        C3 = zeros(2,2,2,2);
        C3(1,:,:,1) = a2;
        C3(2,:,:,2) = a2;
        
        C4 = zeros(2,2,2,2);
        C4(:,1,1,:) = a2;
        C4(:,2,2,:) = a2;
        
        C5 = zeros(2,2,2,2);
        C5(:,1,:,1) = a2;
        C5(:,2,:,2) = a2;
        
        C6 = zeros(2,2,2,2);
        C6(:,:,1,1) = a2;
        C6(:,:,2,2) = a2;
        
        D1 = zeros(2,2,2,2);
        D1(1,1,:,:) = p_outer;
        D1(2,2,:,:) = p_outer;
        
        D2 = zeros(2,2,2,2);
        D2(1,:,1,:) = p_outer;
        D2(2,:,2,:) = p_outer;
        
        D3 = zeros(2,2,2,2);
        D3(1,:,:,1) = p_outer;
        D3(2,:,:,2) = p_outer;
        
        D4 = zeros(2,2,2,2);
        D4(:,1,1,:) = p_outer;
        D4(:,2,2,:) = p_outer;
        
        D5 = zeros(2,2,2,2);
        D5(:,1,:,1) = p_outer;
        D5(:,2,:,2) = p_outer;
        
        D6 = zeros(2,2,2,2);
        D6(:,:,1,1) = p_outer;
        D6(:,:,2,2) = p_outer;
        
        I1 = zeros(2,2,2,2);
        I1(1,1,1,1) = 1;
        I1(2,2,2,2) = 1;
        I1(1,1,2,2) = 1;
        I1(2,2,1,1) = 1;
        
        I2 = zeros(2,2,2,2);
        I2(1,1,1,1) = 1;
        I2(2,2,2,2) = 1;
        I2(1,2,1,2) = 1;
        I2(2,1,1,2) = 1;
        
        I3 = zeros(2,2,2,2);
        I3(1,1,1,1) = 1;
        I3(2,2,2,2) = 1;
        I3(1,2,2,1) = 1;
        I3(2,1,1,2) = 1;
        
        b4 = a4 - (C1 + C2 + C3 + C4 + C5 + C6)/6 + (I1 + I2 + I3)/24;
        
        p4 = zeros(2,2,2,2);
        for i = 1:2
            for j = 1:2
                for k = 1:2
                    for n = 1:2
                        p4(i,j,k,n) = p(i)*p(j)*p(k)*p(n);
                    end
                end
            end
        end
        
        f4 = p4 - (D1 + D2 + D3 + D4 + D5 + D6)/6 +(I1 + I2 + I3)/24;
        
        prob = prob + (8/pi)*sum(b4(:).*f4(:));
    end
          
end %post : recover_probability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = evaluateDLP(o, geom, eta, x, y)
% evaluates the Stokes double layer potential at a point x, y, given a 
% geometry geom and a density function eta
% TO DO : x, y can be a list of points


geomTar = capsules([],[x;y]);

pot = poten(geom.N);
[~,NearStruct] = geom.getZone(geomTar,2);

D = pot.stokesDLmatrix(geom);
DLP = @(X) pot.exactStokesDLdiag(geom,D,X) - 1/2*X;
u = pot.nearSingInt(geom, eta, DLP,[],...
    NearStruct, @pot.exactStokesDL, @pot.exactStokesDL, geomTar,false,false);

end % post : evaluateDLP

end %methods

end %classdef

