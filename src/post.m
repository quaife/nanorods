classdef post
% Post-processing of data from output files

properties
dataFile       

semimajors;
semiminors;
centres_x;
centres_y
orientations;
times;
u_x;
u_y;
omega;
nv;

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = post(dataFile)

o.dataFile = dataFile;

%first 5 lines are header lines
M = dlmread(o.dataFile, '\t', 5, 0);

[~, nc] = size(M);

o.nv = (nc - 1)/6;
o.semimajors = nonzeros(M(1,1:o.nv));
o.semiminors = nonzeros(M(2,1:o.nv));

o.times = M(3:end,1);
o.centres_x = M(3:end,2:o.nv+1);
o.centres_y = M(3:end,o.nv+2:2*o.nv+1);
o.orientations = M(3:end,2*o.nv+2:3*o.nv+1);
o.u_x = M(3:end, 3*o.nv+2:4*o.nv+1);
o.u_y = M(3:end, 4*o.nv+2:5*o.nv+1);
o.omega = M(3:end, 5*o.nv+2:end);


end %post : constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = animated_gif(o, gname, stride, itmax)
    
prams.N = 128;
prams.semimajors = o.semimajors';
prams.semiminors = o.semiminors';

prams.nv = o.nv;

h = figure();

%% find axes limits

if isempty(itmax)
    itmax = length(o.times);
end

xmin = min(min(o.centres_x(1:itmax,:))) - max(max(o.semimajors), max(o.semiminors));
xmax = max(max(o.centres_x(1:itmax,:))) + max(max(o.semimajors), max(o.semiminors));

ymin = min(min(o.centres_y(1:itmax,:))) - max(max(o.semimajors), max(o.semiminors));
ymax = max(max(o.centres_y(1:itmax,:))) + max(max(o.semimajors), max(o.semiminors));


for i = 1:stride:itmax
    
    clf;
    
    xmin = min(min(o.centres_x(i,:))) - max(max(o.semimajors), max(o.semiminors));
    xmax = max(max(o.centres_x(i,:))) + max(max(o.semimajors), max(o.semiminors));
    
%     ymin = min(min(o.centres_y(i,:))) - max(max(o.semimajors), max(o.semiminors));
%     ymax = max(max(o.centres_y(i,:))) + max(max(o.semimajors), max(o.semiminors));
    
    geom = capsules(prams, [o.centres_x(i,:); o.centres_y(i,:)], o.orientations(i,:));
    X = geom.getXY();
    oc = curve;
    [x,y] = oc.getXY(X);
    fill([x;x(1,:)],[y;y(1,:)],'k');
    
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
    axis equal
    
    title(sprintf('t = %6.3f', o.times(i)));
    drawnow
    
    frame = getframe(h);
    im = frame2im(frame);
    
    [imind,cm] = rgb2ind(im,256);
    if i == 1;
        
        imwrite(imind,cm,gname,'gif', 'Loopcount',inf, 'DelayTime', 0);
        
    else
        
        imwrite(imind,cm,gname,'gif','WriteMode','append', 'DelayTime',0);
    end
end

end %post : aimated_gif

function stats = calculate_stats(o, fibers)
%provide statistics for fibers

stats.time = o.times;


stats.semimajors = o.semimajors(fibers);
stats.semiminors = o.semiminors(fibers);
stats.aspect_ratios = o.semiminors(fibers)./o.semimajors(fibers);

stats.centre_x = o.centres_x(:,fibers);
stats.centre_y = o.centres_y(:,fibers);
stats.orientations = wrapToPi(o.orientations(:,fibers));
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

%pdf and tensors
%apply symmetry around 0
p = zeros(size(stats.orientations));
for i = 1:length(stats.time)
    for j = fibers
        p(i,j) = stats.orientations(i,j);
        if (stats.orientations(i,j) < 0)
            p(i,j) = p(i,j) + pi;
        end
    end
end

stats.pdf_thetas = linspace(-pi,pi,length(fibers)/2);
stats.pdf_p = hist(p', stats.pdf_thetas)';
%normalize
dtheta = stats.pdf_thetas(2) - stats.pdf_thetas(1);

for i = 1:length(stats.time)
    stats.pdf_p(i,:) = stats.pdf_p(i,:)/(dtheta*sum(stats.pdf_p(i,:)));
end

%calculate second order moment
stats.a2 = zeros(2,2,length(stats.time));
for i = 1:length(stats.time)
    for j = 1:length(stats.pdf_thetas)
        ptmp = [cos(stats.pdf_thetas(j)); sin(stats.pdf_thetas(j))];
        stats.a2(:,:,i) = stats.a2(:,:,i) + dtheta*(ptmp*ptmp')*stats.pdf_p(i,j);
    end
end

end %post : stats

end %methods

end %classdef

