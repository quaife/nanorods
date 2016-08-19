classdef monitor
% Used for doing input/output of data, monitoring runs during execution

properties
verbose         % write data to console
profile         % profile code (TO DO)
saveData        % save data to the dat files and log file
dataFile        % name of data file containing fibre centres and orientations
densityFile     % name of data file containing density function
logFile         % name of log file
geomFile        % name of geometry file
profileFile     % name of profile folder
append          % append new data to files

OUTPUTPATH_DATA % folder in which to save data
OUTPUTPATH_LOG  % folder in which to save logs
OUTPUTPATH_PROFILE  % folder in which to save logs

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = monitor(options, prams)
% monitor(options,prams) saves options and parameters needed by the class.
% This is the constructor

o.OUTPUTPATH_DATA = '../output/data/';
o.OUTPUTPATH_LOG = '../output/logs/';
o.OUTPUTPATH_PROFILE = '../output/profile/';

o.verbose = options.verbose;
% write data to console

o.saveData = options.saveData;
o.dataFile = [o.OUTPUTPATH_DATA, options.fileBase, '.dat'];
o.densityFile = [o.OUTPUTPATH_DATA, options.fileBase, '_density.dat'];
o.geomFile = [o.OUTPUTPATH_DATA, options.fileBase, '_geometry', '.dat'];
o.logFile = [o.OUTPUTPATH_LOG, options.fileBase, '.log'];
o.profileFile = [o.OUTPUTPATH_PROFILE, options.fileBase];

o.append = options.append;
o.profile = options.profile;

o.welcomeMessage();
%% start new data file if needed
if (o.saveData && ~o.append)
    o.clearFiles();
    fid = fopen(o.dataFile,'w');
    
    fprintf(fid, '%s\n', ['Data file for nanorod simulation ', datestr(now)]);    

    fprintf(fid, '%s\n', 'Line 8 contains length of each rectangular rod');
    fprintf(fid, '%s\n', 'Line 9 contains width of each rectangular rod');
    fprintf(fid, '%s\n', ['Line 10 contains the order that determines the ',...
                        'curvature of the rods']);
    fprintf(fid, '%s\n', ['Line 11 contains: time step order, '...
                         'FMM status, near singular integration status']);
    fprintf(fid, '%s\n', ['Lines 12 onward contain the time, x centre ',...
                        'coordinate, y centre coordinate and the ', ...
                        'orientation for each rod at time t']);

    fprintf(fid,'%s\n', 'BEGIN DATA');
    fprintf(fid,'%s\n', num2str(prams.lengths));
    fprintf(fid,'%s\n', num2str(prams.widths));
    fprintf(fid,'%s\n', num2str(prams.order));
    fprintf(fid,'%s\n', [num2str(options.tstep_order), ' ', ...
                        num2str(double(options.ifmm)), ' ',...
                        num2str(double(options.inear))]);
    fclose(fid);   
end

end % constructor: monitor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clearFiles(o)
% clearFiles() clears the previous log and data files so that there is 
% nothing from previous runs

fid1 = fopen(o.dataFile,'w');
fid2 = fopen(o.densityFile,'w');
fid3 = fopen(o.logFile,'w');
fid4 = fopen(o.geomFile,'w');

fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);

end % clearFiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function welcomeMessage(o)
% welcomeMessage(options,prams) writes specs from the simulation to the
% log file and console

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end

o.writeStars;
message1 = ['RIGID FIBRE SIMULAION ', datestr(now)];
message2 = [' OMP_NUM_THREADS = ', ...
        num2str(getenv('OMP_NUM_THREADS')), ', MATLAB WORKERS = ', num2str(poolsize)];
o.writeMessage(message1);
o.writeMessage(message2);
o.writeStars;

end % welcomeMessage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function restartMessage(o)
% welcomeMessage(options,prams) writes specs from the simulation to the
% log file and console

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end

o.writeStars;
message1 = ['RESTARTING RIGID FIBRE SIMULAION ', datestr(now)];
message2 = [' OMP_NUM_THREADS = ', ...
        num2str(getenv('OMP_NUM_THREADS')), ', MATLAB WORKERS = ', num2str(poolsize)];
o.writeMessage(message1);
o.writeMessage(message2);
o.writeStars;

end % restartMessage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeStars(o)
% writeStars writes a message of stars to the console
% and the log file depending on verbose and saveData

messageStars = '**********************************************************************';
o.writeMessage(messageStars,'%s\n')

end % writeStars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeMessage(o,message,format)
% function writeMessage(message,format) appends message 
% to o.fileName with format

if nargin == 2
  format = '%s\n';
end
% if user doesn't give format, take it to be string followed by a new
% line

if o.saveData
  fid = fopen(o.logFile, 'a');
  fprintf(fid,format,message);
  fclose(fid);
end
% save to log file
if o.verbose
  disp(message)
end
% write to console if verbose==true


end % writeMessage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData(o, prams, t, c, tau, Ux, Uy, omega)
% writeData(t, cp, tau) writes centre point and orientation of each fibre
% to a csv file. This file can be read in to Matlab later for
% postprocessing.

if ~isempty(prams.tracker_fnc)
   p = prams.tracker_fnc(t);
end
fid = fopen(o.dataFile,'a');
fprintf(fid,'%s\n', num2str([t, c(1,:), c(2,:), tau, Ux, Uy, omega, p']));
fclose(fid);

end % writeData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeGeometry(o, walls)
   
    dlmwrite(o.geomFile,walls.X);
end % writeGeometry

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeDensity(o,t, eta)
% writeData(t, eta) write the density function eta at time t to
% o.densityFile. This can be read in later to plot the velocity field.

    fid = fopen(o.densityFile,'a');
    fprintf(fid,'%s\n', num2str([t,eta(:)']));
    fclose(fid);
end % writeDensity

end % methods


end % classdef