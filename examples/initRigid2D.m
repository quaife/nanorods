function [options,prams] = initRigid2D(options,prams)
% set the path and assign options and prams to default values if they
% have not been assigned

% P = path; ii = find(pwd == filesep); ii = ii(end);
% subPath = pwd; subPath = [subPath(1:ii) 'src'];
% if isempty(strfind(P, subPath))
%   addpath(subPath)
% end

PramList = {'Np','np', 'Nw', 'nw', 'lengths', 'widths', 'rounding_order', ...
                        'T', 'number_steps'};
defaultPram.Np = 64;
defaultPram.np = 1;
defaultPram.Nw = 0;
defaultPram.nw = 0;
defaultPram.T = 1;
defaultPram.number_steps = 10;
defaultPram.rounding_order = 2;
for k = 1:length(PramList)
  if ~isfield(prams,PramList{k})
    eval(['prams.' PramList{k} '=defaultPram.' PramList{k} ';'])
    % Set any unassigned parameters to a default value
  end
end


OptionList = {'tstep_order','near_singular','far_field','verbose','use_precond','fmm', ...
        'append', 'gmres_tol', 'profile', 'save_data', 'confined', 'rk_tol', 'rk_min_dt', ...
        'rk_max_up', 'rk_max_down', 'rk_safety'};
    
defaultOption.tstep_order = 2;
defaultOption.near_singular = true;
defaultOption.far_field = 'shear';
defaultOption.verbose = true;
defaultOption.use_precond = false;
defaultOption.fmm = false;
defaultOption.append = false;
defaultOption.profile = false;
defaultOption.gmres_tol = 1e-6;
defaultOption.confined = false;
defaultOption.save_data = true;
defaultOption.rk_tol = 1e-2;
defaultOption.rk_min_dt = 1e-5;
defaultOption.rk_max_up = 1.5;
defaultOption.rk_max_down = 0.5;
defaultOption.rk_safety = 0.9;

for k = 1:length(OptionList)
  if ~isfield(options,OptionList{k})
    eval(['options.' OptionList{k} '=defaultOption.' ...
        OptionList{k} ';'])
    % Set any unassigned parameters to a default value
  end
end


