function setup_code(varargin)
% Get repository path
repo_path = which('setup_dsopt_code.m');
repo_path = fileparts(repo_path);

% remove any auxiliary folder from the search path
restoredefaultpath();

% remove the default user-specific path
userpath('clear');

% Add sub-directories to path
addpath(genpath(repo_path));

% Add libraries to path
libraries_path = strcat(repo_path, '/ds-libraries/');
addpath(genpath(libraries_path));

% Check if compilation is necessary
compile_thirdparty = 1;
if nargin == 1
   compile_thirdparty = varargin{1}; 
end
    
if compile_thirdparty
    
    % NOT WORKING FOR WINDOWS PC (MEX with a compiler needs to be configured = install compiler)
    % Install lightspeed for Chapter 3 and 5 (PC-GMM Fitting)
    lightspeed_path = strcat(repo_path, '/ds-libraries/thirdparty/lightspeed');
    cd(lightspeed_path)
    pause(1);

    % If error here run "install_lightspeed", mex compiler has to be setup "mex -setup"
    install_lightspeed 
    test_lightspeed

    % Go back to ch3 directory
    repo_path = which('setup_ds_code.m');
    repo_path = fileparts(repo_path);
    cd(repo_path)
    
end

% If error here run "install_sedumi -build"
install_sedumi

display('--> Code and toolboxes installed correctly.');
end
