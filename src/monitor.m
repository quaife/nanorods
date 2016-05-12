classdef monitor
% Used for doing input/output of data.
% Can do plotting, write to files, compute error in area
% and length, and write to console

properties
verbose         % write data to console
usePlot         % plot the bodies
%saveData        % save data to the dat files and log file
axis            % axis of the plot
%dataFile        % name of the file to write the data
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
%o.saveData = options.saveData;
%% save messages to a log file
%o.dataFile = options.dataFile;
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


end % constructor: monitor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clearFiles(o)
% clearFiles() clears the log and bin file so that there is nothing from
% previous runs

fid = fopen(o.logFile,'w');
fclose(fid);
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
  fid = fopen(o.logFile,'a');
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
plot([x;x(1,:)],[y;y(1,:)],'k','linewidth',2);
axis equal;
axis(o.axis);

end % plotData


end % methods


end % classdef

