function setup_dsopt_code(varargin)
% Get repository path
repo_path = which('setup_dsopt_code.m');
repo_path = fileparts(repo_path);

% remove any auxiliary folder from the search path
restoredefaultpath();

% remove the default user-specific path
userpath('clear');

% Add sub-directories to path
addpath(genpath(repo_path));

% Add dependent libraries to path (assuming phys-gmm 
libraries_path = strcat(repo_path, '/../phys-gmm/');
addpath(genpath(libraries_path));

% If error here run "install_sedumi -build"
install_sedumi

display('--> Code and toolboxes installed correctly.');
end
