function [options,prams] = initRigid2D(options,prams)
% set the path and assign options and prams to default values if they
% have not been assigned

P = path; ii = find(pwd == filesep); ii = ii(end);
subPath = pwd; subPath = [subPath(1:ii) 'src'];
if isempty(strfind(P, subPath))
  addpath(subPath)
end

PramList = {'N','nv','order', 'T', 'm','gmresTol'};
defaultPram.N = 64;
defaultPram.nv = 1;
defaultPram.T = 1;
defaultPram.m = 10;
defaultPram.gmresTol = 1e-6;
defaultPram.order = 4;
for k = 1:length(PramList)
  if ~isfield(prams,PramList{k})
    eval(['prams.' PramList{k} '=defaultPram.' PramList{k} ';'])
    % Set any unassigned parameters to a default value
  end
end


OptionList = {'tstep_order','inear','farField','verbose','usePreco','ifmm', ...
        'append', 'profile'};
defaultOption.tstep_order = 1;
defaultOption.inear = true;
defaultOption.farField = 'shear';
defaultOption.verbose = true;
defaultOption.usePreco = false;
defaultOption.fmm = false;
defaultOption.append = false;
defaultOption.profile = false;

for k = 1:length(OptionList)
  if ~isfield(options,OptionList{k})
    eval(['options.' OptionList{k} '=defaultOption.' ...
        OptionList{k} ';'])
    % Set any unassigned parameters to a default value
  end
end


