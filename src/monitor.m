classdef monitor
% Used for doing input/output of data, monitoring runs during execution

properties
verbose         % write data to console
profile         % profile code (TO DO)
saveData        % save data to the dat files and log file
dataFile        % name of data file containing fibre centres and orientations
logFile         % name of log file
profileFile     % name of profile folder
append          % append new data to files
options         % options structure
prams           % prameter structure

OUTPUTPATH_DATA % folder in which to save data
OUTPUTPATH_LOG  % folder in which to save logs
OUTPUTPATH_PROFILE  % folder in which to save logs

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = monitor(options, prams, xc, tau)
% monitor(options,prams) saves options and parameters needed by the class.
% This is the constructor

o.OUTPUTPATH_DATA = '../output/data/';
o.OUTPUTPATH_LOG = '../output/logs/';
o.OUTPUTPATH_PROFILE = '../output/profile/';

o.verbose = options.verbose;
% write data to console

o.saveData = options.saveData;
o.dataFile = [o.OUTPUTPATH_DATA, options.fileBase, '.mat'];
o.logFile = [o.OUTPUTPATH_LOG, options.fileBase, '.log'];
o.profileFile = [o.OUTPUTPATH_PROFILE, options.fileBase];

o.append = options.append;
o.profile = options.profile;

o.options = options;
o.prams = prams;

o.welcomeMessage();
%% start new data file if needed
if (o.saveData && ~o.append)

    o.clearFiles();
    etaW = zeros(2*prams.Nbd, prams.nbd);
    etaF = zeros(2*prams.N, prams.nv);
    
    U = zeros(2, prams.nv);
    omega = zeros(1,prams.nv);
    
    stokes = zeros(2,prams.nbd-1);
    rot = zeros(1,prams.nbd-1);
    t = 0;
    
    save(o.dataFile, 'prams', 'options', 't', 'xc', 'tau', 'U', 'omega', 'stokes', 'rot', 'etaW', 'etaF');
end

end % constructor: monitor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clearFiles(o)
% clearFiles() clears the previous log and data files so that there is 
% nothing from previous runs

fid1 = fopen(o.dataFile,'w');
fid2 = fopen(o.logFile,'w');

fclose(fid1);
fclose(fid2);

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
function writeData(o, t_c, xc_c, tau_c, U_c, omega_c, stokes_c, rot_c, etaF_c, etaW_c)
% writeData(t, xc, tau, U, omega, etaF, etaW) writes centre point and orientation of each fibre
% to a mat file. This file can be read in to Matlab later for
% postprocessing.

load(o.dataFile);

t = [t, t_c];

% variables solved at t_(k+1)
xc(:,:,end+1) = xc_c;
tau = [tau; tau_c];


% variables solved at time level t_k
if (length(t) > 2)
    U(:,:,end+1) = U_c;
    omega = [omega;omega_c];
    etaF(:,:,end+1) = etaF_c;

    if o.options.confined
        etaW(:,:,end+1) = etaW_c;

        if o.prams.nbd > 1
            stokes(:,:,end+1) = stokes_c;
            rot(:,end+1) = rot_c;
        end
    end
else
    U(:,:,1) = U_c;
    omega = omega_c;
    etaF(:,:,1) = etaF_c;

    if o.options.confined
        etaW(:,:,1) = etaW_c;

        if o.prams.nbd > 1
            stokes(:,:,1) = stokes_c;
            rot(:,1) = rot_c;
        end
    end
end

save(o.dataFile, 'prams', 'options', 't', 'xc', 'tau', 'U', 'omega', 'stokes', 'rot', 'etaW', 'etaF');

end % writeData

end % methods


end % classdef