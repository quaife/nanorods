clear all
addpath ../src

prams.Nouter = 256;
% number of points on outer solid wall
prams.Ninner = 128;
% number of points per circle exclusion
prams.nv = 3;
% number of exclusions
prams.gmresTol = 1e-6;
% gmres tolerance
prams.maxIter = 100;
% maximum number of gmres iterations

% Different options
options.dataFile = 'output/example1.bin';
options.farField = 'pipe';
options.logFile = 'output/example1.log';
options.saveData = true;

oc = curve;
Xouter = oc.initConfig(prams.Nouter,'square');
% outer most boundary
center = [[-0.3 -0.4];[0.4 0.5];[-0.2 0.2]];
orientation = [pi/3;5*pi/4;2*pi/3];
Xinner = oc.initConfig(prams.Ninner,'rod', ...
          'nv',prams.nv, ...
          'center',center, ...
          'orientation',orientation, ...
          'aspect_ratio',2);
% inner rods

clf; hold on;
plot(Xouter(1:end/2),Xouter(end/2+1:end),'k')
plot(Xinner(1:end/2,:),Xinner(end/2+1:end,:),'r')

%stokesSolver(Xinner,Xouter,options,prams);
% solve density function and write to .bin files



