%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Demo for testing non-linear DS estimation via diff. matching from Perrin %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 1 - OPTION 1 (DATA LOADING): Load CORL-paper Datasets %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc
%%%%%%%%%%%%%%%%%%%%%%%%% Select a Dataset %%%%%%%%%%%%%%%%%%%%%%%%%%
% 1:     Messy Snake Dataset   (2D)
% 2:     L-shape Dataset       (2D)
% 3:     A-shape Dataset       (2D)
% 4:     S-shape Dataset       (2D)
% 5:     Dual-behavior Dataset (2D)
% 6:     Via-point Dataset     (3D) -- 15 trajectories recorded at 100Hz
% 7:     Sink Dataset          (3D) -- 21 trajectories recorded at 100Hz
% 8:     CShape top            (3D) -- 10 trajectories recorded at 100Hz
% 9:     CShape bottom         (3D) -- 10 trajectories recorded at 100Hz
% 10:    CShape all            (3D) -- 20 trajectories recorded at 100Hz
% 11:    Cube arranging        (3D) -- 20 trajectories recorded at 100Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg_dir         = '/home/nbfigueroa/Dropbox/PhD_papers/CoRL-2018/code/ds-opt/';
chosen_dataset  = 9; 
sub_sample      = 2; % '>2' for real 3D Datasets, '1' for 2D toy datasets
nb_trajectories = 10; % For real 3D data
[Data, Data_sh, att, x0_all, data, dt] = load_dataset_DS(pkg_dir, chosen_dataset, sub_sample, nb_trajectories);

% Position/Velocity Trajectories
vel_samples = 10; vel_size = 0.5; 
[h_data, h_att, h_vel] = plot_reference_trajectories_DS(Data, att, vel_samples, vel_size);
dim =size(data{1},1)/2;
nb_demos    = length(data);

%% %%%%%%%%%%%% [Optional] Load pre-learned SEDS model from Mat file  %%%%%%%%%%%%%%%%%%%
DS_name = '3D-Sink_diff';
matfile = strcat(pkg_dir,'/models/', DS_name,'.mat');
load(matfile)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 1 - OPTION 2 (DATA LOADING): Load Motions from LASA Handwriting Dataset %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose DS LASA Dataset to load
clear all; close all; clc

% Select one of the motions from the LASA Handwriting Dataset
sub_sample      = 1; % Each trajectory has 1000 samples when set to '1'
nb_trajectories = 5; % Maximum 7, will select randomly if <7
[Data, Data_sh, att, x0_all, data, dt] = load_LASA_dataset_DS(sub_sample, nb_trajectories);

% Position/Velocity Trajectories
vel_samples = 15; vel_size = 0.5; 
[h_data, h_att, h_vel] = plot_reference_trajectories_DS(Data, att, vel_samples, vel_size);
dim =size(data{1},1)/2;
nb_demos    = length(data);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2: Estimate Diffeomorphic Matching Function Parameters  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some options
step_demos  = 1;  %Step the points
doNormalize = 0;  % Normalize the demos using the variance in each direction; This works better in general
                  % Nadia's comment: This indeed works better but deforms
                  % the demonstrations!

%Options concerning the transformation (NOTE:: Normalization deforms the trajectories!)
if doNormalize
    optsSearch = {'maxCoef', '10', 'nb_iteration','150', 'regularise', '5e-4', 'conv_crit', '1e-6', 'division_coefList', '[3,3,2.5,2.0,1.5,1.1,1.1]', 'safeCoeffList', '0.6*[1, 1, 1, 1, 1, 1, 1]', 'doPlot', '0'};
else
    optsSearch = {'maxCoef', '5', 'nb_iteration','150', 'regularise', '1e-3', 'conv_crit', '1e-6', 'division_coefList', '[3,3,2.5,2.0,1.5,1.1,1.1]', 'safeCoeffList', '0.55*[1, 1, 1, 1, 1, 1, 1]', 'doPlot', '0'};
end

% Prepare data
fprintf('- Preparing Data for Diffeomorphic Matching...');
[Xinit, target_trajectory, alltarget_trajectory, alltarget_trajectoryV, allSource,  ... ,
allSourceV, nDemos, indD] = processDrawnDataset(data, 1:nb_demos, step_demos,dt);
fprintf('done\n');

% Uses the mean of all demonstrations as target_trajectory Xi
% (The transformation maps the points given in "source" onto the points in "target_trajectory"
[dim, lX] = size(target_trajectory);
target_trajectory = target_trajectory./nDemos;
target_trajectoryV = [diff(target_trajectory,[],2),zeros(dim,1)]./dt;

%Define the source trajectory: A straight line between the initial and
%final point of the mean target_trajectory trajectory
if dim == 2
    source = [linspace(target_trajectory(1,1), target_trajectory(1,end), lX);
              linspace(target_trajectory(2,1), target_trajectory(2,end), lX)];
elseif dim == 3
    source = [linspace(target_trajectory(1,1), target_trajectory(1,end), lX);
              linspace(target_trajectory(2,1), target_trajectory(2,end), lX);
              linspace(target_trajectory(3,1), target_trajectory(3,end), lX);];
end
sourceV = [diff(source, [], 2).*(step_demos/dt), zeros(dim,1)];

% Sub-sample source/target_trajectory trajectories depending on the step-size
source  = source(:,indD);
sourceV = sourceV(:,indD);
target_trajectory  = target_trajectory(:,indD);
target_trajectoryV = target_trajectoryV(:,indD);

% Normalize Trajectories
if doNormalize
    [Xinit, target_trajectory, target_trajectoryV, source, sourceV, alltarget_trajectory, ... ,
        alltarget_trajectoryV, allSource, allSourceV] = normalizeTrajectories(Xinit, target_trajectory, ... ,
        target_trajectoryV, source, sourceV, alltarget_trajectory, alltarget_trajectoryV, allSource, allSourceV);
end

%Search for the transformation parameters
fprintf('- Searching for Transformation Parameters...');
tic;
[ centers, target_trajectorys, coefs, division_coef, nb_iteration ] = iterativeSearch( source, target_trajectory, optsSearch{:} );
toc;
fprintf('done\n')

% Create all lambda functions needed to evaluate the transformation
% Forward transformation (Phi)
diff_fun             = @(pt) result_function(centers, target_trajectorys, coefs, division_coef, nb_iteration, pt);
% Backward transformation (Phi^-1)
inverse_diff_fun     = @(pt) result_function_reverse(centers, target_trajectorys, coefs, division_coef, nb_iteration, pt);
% Backward transformation and calculation of the jacobian in point (J^-1)
jac_inverse_diff_fun = @(pt) result_function_reverse_Jac(centers, target_trajectorys, coefs, division_coef, nb_iteration, pt);

%%%%%%%%%%%% General plot with the transformed sources and target_trajectorys %%%%%%%%%%%%
close all;
h_fig = plotDiffeomorphicTransformation(nDemos, source, target_trajectory, alltarget_trajectory, allSource, diff_fun, inverse_diff_fun);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Step 3: Generated Deformed Dynamics function       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate DS function
lambda_1 = 1;
EIG0 = -diag([lambda_1,2*lambda_1.*ones(1,dim-1)]);
ds_diff = @(x) diffeomorphic_ds(x-repmat(att,[1 size(x,2)]), EIG0, source, jac_inverse_diff_fun);

%%%%%%%%%%%%%%    Plot Resulting DS  %%%%%%%%%%%%%%%%%%%
% Fill in plotting options
ds_plot_options = [];
ds_plot_options.sim_traj  = 1;            % To simulate trajectories from x0_all
ds_plot_options.x0_all    = x0_all;       % Intial Points
ds_plot_options.init_type = 'cube';       % For 3D DS, to initialize streamlines
                                          % 'ellipsoid' or 'cube'  
ds_plot_options.nb_points = 30;           % No of streamlines to plot (3D)
ds_plot_options.plot_vol  = 1;            % Plot volume of initial points (3D)
[hd, hs, hr, x_sim] = visualizeEstimatedDS(Data(1:dim,:), ds_diff, ds_plot_options);
limits = axis;
title('Diffeomorphic Matching - Dynamics', 'Interpreter','LaTex','FontSize',17)

%% %%%%%%%%%%%%   Export diff-DS parameters to Mat file  %%%%%%%%%%%%%%%%%%%
DS_name = '3D-CShape-top_diff';
save_diffDS_to_Mat(DS_name, pkg_dir, EIG0, source, jac_inverse_diff_fun,  att, x0_all, dt)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 4 (Evaluation): Compute Metrics and Visualize Velocities %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% NOTE: THIS WILL NOT WORK IF YOU NORMALIZE
%% Compute Errors
% Compute RMSE on training data
rmse = mean(rmse_error(ds_diff, Data(1:dim,:), Data(dim+1:end,:)));
fprintf('diff-DS got prediction RMSE on training set: %d \n', rmse);

% Compute e_dot on training data
edot = mean(edot_error(ds_diff, Data(1:dim,:), Data(dim+1:end,:)));
fprintf('diff-DS got prediction e_dot on training set: %d \n', edot);

% Compute DTWD between train trajectories and reproductions
if ds_plot_options.sim_traj
    nb_traj       = size(x_sim,3);
    ref_traj_leng = size(Data,2)/nb_traj;
    dtwd = zeros(1,nb_traj);
    for n=1:nb_traj
        start_id = round(1+(n-1)*ref_traj_leng);
        end_id   = round(n*ref_traj_leng);
        dtwd(1,n) = dtw(x_sim(:,:,n)',Data(1:dim,start_id:end_id)',20);
    end
    fprintf('diff-DS got DTWD of reproduced trajectories: %2.4f +/- %2.4f \n', mean(dtwd),std(dtwd));
end
% Compare Velocities from Demonstration vs DS
h_vel = visualizeEstimatedVelocities(Data, ds_diff);
