classdef monitor
% Used for doing input/output of data.
% Can do plotting, write to files, compute error in area
% and length, and write to console

properties
verbose         % write data to console
usePlot         % plot the bodies
saveData        % save data to the dat files and log file
axis            % axis of the plot
dataFile        % name of the file to write the data
append
%logFile         % name of the file to write the log
N               % number of points on inner boundaries
%Nouter          % number of points on outer boundary
%Ntracers        % number of tracers
nv              % number of inner boundaries

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = monitor(options,prams)
% monitor(X,Xwalls,options,prams) saves options and parameters
% needed by the class and also computes the initial error in
% area and length so that the errors in area and length
% can be monitored throughout the simulation.
% This is the constructor


o.verbose = options.verbose;
% write data to console
o.usePlot = options.usePlot;
% plot the rigid bodies
o.saveData = options.saveData;
%% save messages to a log file
o.dataFile = options.dataFile;
o.append = options.append;
%% name of bin file for geometry
%% Uses this file name to choose the file name for the tracers
%o.logFile = options.logFile;
%% name of log file
o.axis = options.axis;
% axis for the plotting
o.N = prams.N;
% number of points per inner boundary
o.nv = prams.nv;
% number of inner boundaries

%% start new data file if needed
if (o.saveData && ~o.append)
    o.clearFiles();
    
    o.writeMessage('Data file for nanorod simulation');
    o.writeMessage('Line 6 contains semi-major axes for each elliptical rod')
    o.writeMessage('Line 7 contains semi-minor axes for each elliptical rod');
    o.writeMessage('Lines 8 onward contain the time, x centre coordinate, y centre coordinate and the orientation for each rod at time t');
    o.writeMessage('BEGIN DATA');
    
    fid = fopen(o.dataFile,'a');
    fprintf(fid,'%s\n', num2str(prams.semimajor));
    fprintf(fid,'%s\n', num2str(prams.semiminor));
    fclose(fid);
end

end % constructor: monitor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clearFiles(o)
% clearFiles() clears the log and bin file so that there is nothing from
% previous runs

% fid = fopen(o.logFile,'w');
% fclose(fid);
fid = fopen(o.dataFile,'w');
fclose(fid);
% delete the previous log and data files

end % clearFiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function welcomeMessage(o,options,prams)
% welcomeMessage(options,prams) writes specs from the simulation to the
% log file and console

o.writeStars
message = ['stuff '];
o.writeMessage(message);

end % welcomeMessage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeStars(o)
% writeStars writes a message of stars to the console
% and the log file depending on verbose and saveData

messageStars = '*********************************************';
o.writeMessage(messageStars,'%s\n')

end % writeStars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeMessage(o,message,format)
% function writeMessage(message,format) appends message 
% to o.fileName with format

if nargin == 2
  format = '%s\n';
end
% if user doesn't give format, take it to be string followed by a new
% line

if o.saveData
  fid = fopen(o.dataFile,'a');
  fprintf(fid,format,message);
  fclose(fid);
end
% save to log file
if o.verbose
  disp(message)
end
% write to console if verbose==true


end % writeMessage


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outputInfo(o,X)

if o.usePlot
  o.plotData(X); 
  pause(1e-2);
end

end % outputInfo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotData(o,X);

oc = curve;
[x,y] = oc.getXY(X);
figure(1);clf;
fill([x;x(1,:)],[y;y(1,:)],'k');
axis equal;
axis(o.axis);

end % plotData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData(o, t, c, tau)
% writeData(t, cp, tau) writes centre point and orientation of each fibre
% to a csv file. This file can be read in to Matlab later for
% postprocessing.

fid = fopen(o.dataFile,'a');
fprintf(fid,'%s\n', num2str([t, c(1,:), c(2,:), tau]));
fclose(fid);

end

end % methods


end % classdef

