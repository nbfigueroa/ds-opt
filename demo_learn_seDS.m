%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo Code for SEDS Learning Comparison for paper:                       %
%  'A Physically Const....'                                               %
% Author: Nadia Figueroa                                                  %
% Date: June 3rd, 2018                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 1 - OPTION 1 (DATA LOADING): Load CORL-paper Datasets %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc
%%%%%%%%%%%%%%%%%%%%%%%%% Select a Dataset %%%%%%%%%%%%%%%%%%%%%%%%%%
% 1:  Messy Snake Dataset   (2D)
% 2:  L-shape Dataset       (2D)
% 3:  A-shape Dataset       (2D)
% 4:  S-shape Dataset       (2D)
% 5:  Dual-behavior Dataset (2D)
% 6:  Via-point Dataset     (3D) -- 15 trajectories recorded at 100Hz
% 7:  Sink Dataset          (3D) -- 21 trajectories recorded at 100Hz
% 8:  CShape top            (3D) -- 10 trajectories recorded at 100Hz
% 9:  CShape bottom         (3D) -- 10 trajectories recorded at 100Hz
% 10: CShape all            (3D) -- 20 trajectories recorded at 100Hz
% 11: Cube arranging        (3D) -- 20 trajectories recorded at 100Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg_dir         = '/home/nbfigueroa/Dropbox/PhD_papers/CoRL-2018/code/ds-opt/';
chosen_dataset  = 1; 
sub_sample      = 1; % '>2' for real 3D Datasets, '1' for 2D toy datasets
nb_trajectories = 0; % Only for real 3D data
[Data, Data_sh, att, x0_all, data, ~] = load_dataset_DS(pkg_dir, chosen_dataset, sub_sample, nb_trajectories);

% Position/Velocity Trajectories
vel_samples = 10; vel_size = 0.5; 
[h_data, h_att, h_vel] = plot_reference_trajectories_DS(Data, att, vel_samples, vel_size);
limits = axis;

% Extract Position and Velocities
M          = size(Data,1)/2;    
Xi_ref     = Data_sh(1:M,:);
Xi_dot_ref = Data_sh(M+1:end,:);     

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 1 - OPTION 2 (DATA LOADING): Load Motions from LASA Handwriting Dataset %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose DS LASA Dataset to load
clear all; close all; clc

% Select one of the motions from the LASA Handwriting Dataset
sub_sample      = 3; % Each trajectory has 1000 samples when set to '1'
nb_trajectories = 7; % Maximum 7, will select randomly if <7
[Data, Data_sh, att, x0_all, ~, dt] = load_LASA_dataset_DS(sub_sample, nb_trajectories);

% Position/Velocity Trajectories
vel_samples = 15; vel_size = 0.5; 
[h_data, h_att, h_vel] = plot_reference_trajectories_DS(Data, att, vel_samples, vel_size);

% Extract Position and Velocities
M          = size(Data,1)/2;    
Xi_ref     = Data(1:M,:);
Xi_dot_ref = Data(M+1:end,:);  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 2 (GMM FITTING): Fit GMM to Trajectory Data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'seds-init': follows the initialization given in the SEDS code
% 0: Manually set the # of Gaussians
% 1: Do Model Selection with BIC
do_ms_bic = 1;

if do_ms_bic
    est_options = [];
    est_options.type        = 1;   % GMM Estimation Alorithm Type
    est_options.maxK        = 15;  % Maximum Gaussians for Type 1/2
    est_options.do_plots    = 1;   % Plot Estimation Statistics
    est_options.fixed_K     = [];   % Fix K and estimate with EM
    est_options.sub_sample  = 1;   % Size of sub-sampling of trajectories 
    
    [Priors0, Mu0, Sigma0] = fit_gmm([Xi_ref; Xi_dot_ref], [], est_options);
    nb_gaussians = length(Priors0);
else
    % Select manually the number of Gaussian components
    nb_gaussians = 5;
end

%finding an initial guess for GMM's parameter
[Priors0, Mu0, Sigma0] = initialize_SEDS([Xi_ref; Xi_dot_ref],nb_gaussians);
title_string = '$\theta_{\gamma}=\{\pi_k,\mu^k,\Sigma^k\}$ Initial Estimate';

% Plot Initial Estimate 
[~, est_labels] =  my_gmm_cluster([Xi_ref; Xi_dot_ref], Priors0, Mu0, Sigma0, 'hard', []);
% Position/Velocity Trajectories
vel_samples = 20; vel_size = 0.5; 
[h_data, h_att, h_vel] = plot_reference_trajectories_DS(Data_sh, [0 0]', vel_samples, vel_size);
limits = axis;
plotGMMParameters( Xi_ref, est_labels, Mu0(1:2,:), Sigma0(1:2,1:2,:),1);
title(title_string, 'Interpreter', 'LaTex', 'FontSize',20)
ml_plot_gmm_pdf(Xi_ref, Priors0, Mu0(1:2,:), Sigma0(1:2,1:2,:), limits)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%  Step 3 (DS ESTIMATION): RUN SEDS SOLVER  %%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear options;
options.tol_mat_bias  = 10^-6;    % A very small positive scalar to avoid
                                  % instabilities in Gaussian kernel [default: 10^-1]                             
options.display       = 1;        % An option to control whether the algorithm
                                  % displays the output of each iterations [default: true]                            
options.tol_stopping  = 10^-9;    % A small positive scalar defining the stoppping
                                  % tolerance for the optimization solver [default: 10^-10]
options.max_iter      = 1000;     % Maximum number of iteration forthe solver [default: i_max=1000]
options.objective     = 'mse';    % 'mse'/'likelihood'
sub_sample            = 1;

%running SEDS optimization solver
[Priors Mu Sigma]= SEDS_Solver(Priors0,Mu0,Sigma0,[Xi_ref(:,1:sub_sample:end); Xi_dot_ref(:,1:sub_sample:end)],options); 
ds_seds = @(x) GMR_SEDS(Priors,Mu,Sigma,x-repmat(att,[1 size(x,2)]),1:2,3:4);

%%%%%%%%%%%%%%    Plot Resulting DS  %%%%%%%%%%%%%%%%%%%
simulate_reproductions = 1;
[hd, hs, hr, x_sim] = visualizeEstimatedDS(Data(1:2,:), ds_seds, simulate_reproductions, x0_all);
limits = axis;
switch options.objective
    case 'mse'        
        title('SEDS Dynamics with $J(\theta_{\gamma})$=MSE', 'Interpreter','LaTex','FontSize',20)
    case 'likelihood'
        title('SEDS Dynamics with $J(\theta_{\gamma})$= log-Likelihood', 'Interpreter','LaTex','FontSize',20)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 4 (Evaluation): Compute Metrics and Visualize Velocities %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Errors
clc;
% Compute RMSE on training data
rmse = mean(rmse_error(ds_seds, Xi_ref, Xi_dot_ref));
fprintf('SEDS got prediction RMSE on training set: %d \n', rmse);

% Compute e_dot on training data
edot = mean(edot_error(ds_seds, Xi_ref, Xi_dot_ref));
fprintf('SEDS got prediction e_dot on training set: %d \n', edot);

% Compute DTWD between train trajectories and reproductions
if simulate_reproductions
    nb_traj       = size(x_sim,3);
    ref_traj_leng = size(Xi_ref,2)/nb_traj;
    dtwd = zeros(1,nb_traj);
    for n=1:nb_traj
        start_id = round(1+(n-1)*ref_traj_leng);
        end_id   = round(n*ref_traj_leng);
        dtwd(1,n) = dtw(x_sim(:,:,n)',Data(1:2,start_id:end_id)',20);
    end
end
fprintf('SEDS got reproduction DTWD on training set: %2.4f +/- %2.4f \n', mean(dtwd),std(dtwd));

% Compare Velocities from Demonstration vs DS
h_vel = visualizeEstimatedVelocities(Data, ds_seds);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Step 5 (Optional - Stability Check 2D-only): Plot Lyapunov Function and derivative  %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Type of plot
contour = 1; % 0: surf, 1: contour
clear lyap_fun_comb lyap_der 
P = eye(2);
% Lyapunov function
lyap_fun = @(x)lyapunov_function_PQLF(x, att, P);
title_string = {'$V(\xi) = (\xi-\xi^*)^T(\xi-\xi^*)$'};

% Derivative of Lyapunov function (gradV*f(x))
lyap_der = @(x)lyapunov_derivative_PQLF(x, att, P, ds_seds);
title_string_der = {'Lyapunov Function Derivative $\dot{V}(\xi)$'};

% Plots
h_lyap     = plot_lyap_fct(lyap_fun, contour, limits,  title_string, 0);
[hd] = scatter(Data(1,:),Data(2,:),10,[1 1 0],'filled'); hold on;
h_lyap_der = plot_lyap_fct(lyap_der, contour, limits,  title_string_der, 1);
[hd] = scatter(Data(1,:),Data(2,:),10,[1 1 0],'filled'); hold on;
