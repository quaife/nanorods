classdef post
% Post-processing of data from output files

properties
dataFile       

semimajors;
semiminors;

widths;
lengths;
order;

centres_x;
centres_y
orientations;
times;
u_x;
u_y;
omega;
nv;
capsule_type;

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = post(dataFile, type)

o.capsule_type = type;    
o.dataFile = dataFile;

%first 5 lines are header lines



switch o.capsule_type
    case 'rectangle'

    M = dlmread(o.dataFile, '', 6, 0);

    [~, nc] = size(M);
    o.nv = (nc - 1)/6;
    
    o.lengths = nonzeros(M(1,1:o.nv));
    o.widths = nonzeros(M(2,1:o.nv));
    o.order = nonzeros(M(3,1:o.nv));

    o.times = M(4:end,1);
    o.centres_x = M(4:end,2:o.nv+1);
    o.centres_y = M(4:end,o.nv+2:2*o.nv+1);
    o.orientations = M(4:end,2*o.nv+2:3*o.nv+1);
    o.u_x = M(4:end, 3*o.nv+2:4*o.nv+1);
    o.u_y = M(4:end, 4*o.nv+2:5*o.nv+1);
    o.omega = M(4:end, 5*o.nv+2:end);   
    
    case 'ellipsoid'
        
    M = dlmread(o.dataFile, '', 5, 0);

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
end



end %post : constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = animated_gif(o, gname, stride, itmax)
    
prams.N = 128;
prams.nv = o.nv;
prams.capsule_type = o.capsule_type;
prams.order = o.order;
prams.lengths = o.lengths;
prams.widths = o.widths;

h = figure();

%% find axes limits

if isempty(itmax)
    itmax = length(o.times);
end

switch (o.capsule_type)
    case 'rectangle'
        prams.length = o.lengths;
        prams.width = o.widths;
        prams.order = o.order;
        
        xmin = min(min(o.centres_x(1:itmax,:))) - max(o.lengths);
        xmax = max(max(o.centres_x(1:itmax,:))) + max(o.lengths);
        
        ymin = min(min(o.centres_y(1:itmax,:))) - max(o.lengths);
        ymax = max(max(o.centres_y(1:itmax,:))) + max(o.lengths);
        
    case 'ellipsoid'
        prams.semimajors = o.semimajors';
        prams.semiminors = o.semiminors';
        xmin = min(min(o.centres_x(1:itmax,:))) - max(max(o.semimajors), max(o.semiminors));
        xmax = max(max(o.centres_x(1:itmax,:))) + max(max(o.semimajors), max(o.semiminors));
        
        ymin = min(min(o.centres_y(1:itmax,:))) - max(max(o.semimajors), max(o.semiminors));
        ymax = max(max(o.centres_y(1:itmax,:))) + max(max(o.semimajors), max(o.semiminors));
        
end
for i = 1:stride:itmax
    
    clf;
    
%     switch (o.capsule_type)
%     case 'rectangle'
%         prams.length = o.lengths;
%         prams.width = o.widths;
%         prams.order = o.order;
%         
%         xmin = min(min(o.centres_x(i,:))) - max(o.lengths);
%         xmax = max(max(o.centres_x(i,:))) + max(o.lengths);
%         
%         ymin = min(min(o.centres_y(1:itmax,:))) - max(o.lengths);
%         ymax = max(max(o.centres_y(1:itmax,:))) + max(o.lengths);
%         
%     case 'ellipsoid'
%         prams.semimajors = o.semimajors';
%         prams.semiminors = o.semiminors';
%         xmin = min(min(o.centres_x(1:itmax,:))) - max(max(o.semimajors), max(o.semiminors));
%         xmax = max(max(o.centres_x(1:itmax,:))) + max(max(o.semimajors), max(o.semiminors));
%         
%         ymin = min(min(o.centres_y(1:itmax,:))) - max(max(o.semimajors), max(o.semiminors));
%         ymax = max(max(o.centres_y(1:itmax,:))) + max(max(o.semimajors), max(o.semiminors));
%         
%     end
    
    
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
                        stats.a4(iT,jT,kT,nT,i) = stats.a4(iT,jT,kT,nT,i) + dtheta*ptmp(iT)*ptmp(jT)*ptmp(kT)*ptmp(nT);
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
end %methods

end %classdef

