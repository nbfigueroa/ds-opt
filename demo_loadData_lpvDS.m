%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo Code for GMM-based LPV-DS Learning introduced in paper:            %
%  'A Physically-Consistent Bayesian Non-Parametric Mixture Model for     %
%   Dynamical System Learning.'                                           %
% With this script you can load 2D toy trajectories or even real-world 
% trajectories acquired via kinesthetic taching and test the different    %
% GMM fitting approaches.                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 1 (DATA LOADING): Load Datasets %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
chosen_dataset  = 4; 
sub_sample      = 1; % '>2' for real 3D Datasets, '1' for 2D toy datasets
nb_trajectories = 0; % For real 3D data
[Data, Data_sh, att, x0_all, ~, ~] = load_dataset_DS(pkg_dir, chosen_dataset, sub_sample, nb_trajectories);

% Position/Velocity Trajectories
vel_samples = 10; vel_size = 0.5; 
[h_data, h_att, h_vel] = plot_reference_trajectories_DS(Data, att, vel_samples, vel_size);

% Extract Position and Velocities
M          = size(Data,1)/2;    
Xi_ref     = Data(1:M,:);
Xi_dot_ref = Data(M+1:end,:);       

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 2 (GMM FITTING): Fit GMM to Trajectory Data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% GMM Estimation Algorithm %%%%%%%%%%%%%%%%%%%%%%
% 0: Physically-Consistent Non-Parametric (Collapsed Gibbs Sampler)
% 1: GMM-EM Model Selection via BIC
% 2: CRP-GMM (Collapsed Gibbs Sampler)
est_options = [];
est_options.type             = 1;   % GMM Estimation Alorithm Type   

% If algo 1 selected:
est_options.maxK             = 15;  % Maximum Gaussians for Type 1
est_options.fixed_K          = [];  % Fix K and estimate with EM for Type 1

% If algo 0 or 2 selected:
est_options.samplerIter      = 20;  % Maximum Sampler Iterations
                                    % For type 0: 20-50 iter is sufficient
                                    % For type 2: >100 iter are needed
                                    
est_options.do_plots         = 1;   % Plot Estimation Statistics
est_options.sub_sample       = 1;   % Size of sub-sampling of trajectories

% Metric Hyper-parameters
est_options.estimate_l       = 1;   % '0/1' Estimate the lengthscale, if set to 1
est_options.l_sensitivity    = 1;   % lengthscale sensitivity [1-10->>100]
                                    % Default value is set to '2' as in the
                                    % paper, for very messy, close to
                                    % self-interescting trajectories, we
                                    % recommend a higher value
est_options.length_scale     = [];  % if estimate_l=0 you can define your own
                                    % l, when setting l=0 only
                                    % directionality is taken into account

% Fit GMM to Trajectory Data
[Priors, Mu, Sigma] = fit_gmm(Xi_ref, Xi_dot_ref, est_options);

% Generate GMM data structure for DS learning
clear ds_gmm; ds_gmm.Mu = Mu; ds_gmm.Sigma = Sigma; 
ds_gmm.Priors = Priors; 

%% (Optional) Step 2.1: Adjust the Covariance matrices if they are too thin
% This is particularly useful for EM-estimates
adjusts_C  = 1;
if adjusts_C  == 1
    tot_scale_fact = 1.15; rel_scale_fact = 0.15;
    Sigma_ = adjust_Covariances(Sigma, tot_scale_fact, rel_scale_fact);
    ds_gmm.Sigma = Sigma_;
end   

%%  (Optional) Step 2.2: Visualize Gaussian Components and labels on clustered trajectories 
% Extract Cluster Labels
est_K      = length(Priors);
[~, est_labels] =  my_gmm_cluster(Xi_ref, ds_gmm.Priors, ds_gmm.Mu, ds_gmm.Sigma, 'hard', []);

% Visualize Estimated Parameters
[h_gmm]  = visualizeEstimatedGMM(Xi_ref,  ds_gmm.Priors, ds_gmm.Mu, ds_gmm.Sigma, est_labels, est_options);
limits = axis;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%  Step 3 (DS ESTIMATION): ESTIMATE SYSTEM DYNAMICS MATRICES  %%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% DS OPTIMIZATION OPTIONS %%%%%%%%%%%%%%%%%%%%
% Type of constraints/optimization 
constr_type = 1;      % 0:'convex':     A' + A < 0
                      % 1:'non-convex': A'P + PA < 0
                      % 2:'non-convex': A'P + PA < -Q given P                                  
init_cvx    = 0;      % 0/1: initialize non-cvx problem with cvx                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if constr_type == 0 || constr_type == 1
    P_opt = [];
else
    [Vxf] = learn_wsaqf(Data, att);
    P_opt = Vxf.P(:,:,1);
end

%%%%%%%%  LPV system sum_{k=1}^{K}\gamma_k(xi)(A_kxi + b_k) %%%%%%%%  
if constr_type == 1
    [A_g, b_g, P_est] = optimize_lpv_ds_from_data(Data_sh, [0 0]', constr_type, ds_gmm, P_opt, init_cvx);
    ds_lpv = @(x) lpv_ds(x-repmat(att,[1 size(x,2)]), ds_gmm, A_g, b_g);
else
    [A_g, b_g, P_est] = optimize_lpv_ds_from_data(Data, att, constr_type, ds_gmm, P_opt, init_cvx);
    ds_lpv = @(x) lpv_ds(x, ds_gmm, A_g, b_g);
end

%%%%%%%%%%%%%%    Plot Resulting DS  %%%%%%%%%%%%%%%%%%%
simulate_reproductions = 1;
[hd, hs, hr, x_sim] = visualizeEstimatedDS(Xi_ref, ds_lpv, simulate_reproductions, x0_all);
switch constr_type
    case 0
        title('GMM-based LPV-DS with QLF (O1)', 'Interpreter','LaTex','FontSize',20)
    case 1
        title('GMM-based LPV-DS with P-QLF (O2) ', 'Interpreter','LaTex','FontSize',20)
    case 2
        title('GMM-based LPV-DS with P-QLF (O3)', 'Interpreter','LaTex','FontSize',20)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 4 (Evaluation): Compute Metric and Visualize Velocities  %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Errors
% Compute RMSE on training data
rmse = mean(rmse_error(ds_lpv, Xi_ref, Xi_dot_ref));
fprintf('LPV-DS with (O%d), got prediction RMSE on training set: %d \n', constr_type+1, rmse);

% Compute e_dot on training data
edot = mean(edot_error(ds_lpv, Xi_ref, Xi_dot_ref));
fprintf('LPV-DS with (O%d), got e_dot on training set: %d \n', constr_type+1, edot);

% Compute DTWD between train trajectories and reproductions
if simulate_reproductions
    nb_traj       = size(x_sim,3);
    ref_traj_leng = size(Xi_ref,2)/nb_traj;
    dtwd = zeros(1,nb_traj);
    for n=1:nb_traj
        start_id = round(1+(n-1)*ref_traj_leng);
        end_id   = round(n*ref_traj_leng);
        dtwd(1,n) = dtw(x_sim(:,:,n)',Xi_ref(:,start_id:end_id)',20);
    end
    fprintf('LPV-DS got DTWD of reproduced trajectories: %2.4f +/- %2.4f \n', mean(dtwd),std(dtwd));
end
% Compare Velocities from Demonstration vs DS
h_vel = visualizeEstimatedVelocities(Data, ds_lpv);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Step 5 (Optional - Stability Check): Plot Lyapunov Function and derivative  %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
