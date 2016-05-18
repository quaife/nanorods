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
nv;

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = post(dataFile)

o.dataFile = dataFile;

%first 5 lines are header lines
M = dlmread(o.dataFile, '\t', 5, 0);

[~, nc] = size(M);

o.nv = (nc - 1)/3;
o.semimajors = nonzeros(M(1,1:o.nv));
o.semiminors = nonzeros(M(2,1:o.nv));

o.times = M(3:end,1);
o.centres_x = M(3:end,2:o.nv+1);
o.centres_y = M(3:end,o.nv+2:2*o.nv+1);
o.orientations = M(3:end,2*o.nv+2:end);


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

end

end %methods

end %classdef

