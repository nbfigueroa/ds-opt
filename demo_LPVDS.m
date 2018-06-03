%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo Code for LPV-DS Learning for paper:                                %
%  'A Physically Const....'                                               %
% Author: Nadia Figueroa                                                  %
% Date: June 3rd, 2018                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA LOADING OPTION 1: Draw with GUI %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc
fig1 = figure('Color',[1 1 1]);
% Axis limits
limits = [-6 0.5 -0.5 2];
%     limits = [-6 4 -2 2];
axis(limits)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.55, 0.2646 0.4358]);
grid on

% Global Attractor of DS
att_g = [0 0]';
radius_fun = @(x)(1 - my_exp_loc_act(5, att_g, x));
scatter(att_g(1),att_g(2),100,[0 0 0],'d'); hold on;

% Draw Reference Trajectories
data = draw_mouse_data_on_DS(fig1, limits);
Data = [];
for l=1:length(data)    
    % Check where demos end and shift
    data_ = data{l};
    if radius_fun(data_(1:2,end)) > 0.75
        data_(1:2,:) = data_(1:2,:) - repmat(data_(1:2,end), [1 length(data_)]);
        data_(3:4,end) = zeros(2,1);
    end    
    Data = [Data data_];
end

% Position/Velocity Trajectories
Xi_ref     = Data(1:2,:);
Xi_dot_ref = Data(3:end,:);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA LOADING OPTION 2: Choose from LASA DATASET %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose DS LASA Dataset to load
[demos, limits] = load_LASA_dataset();

% Global Attractor of DS
radius_fun = @(x)(1 - my_exp_loc_act(5, att_g, x));
att_g = [0 0]';

sample = 2;
Data = [];
for l=1:2   
    % Check where demos end and shift    
    data_pos_raw = demos{l}.pos(:,1:sample:end); 
    data_filt    = sgolay_time_derivatives(data_pos_raw', demos{l}.dt, 2, 3, 15);
    data_ = [data_filt(:,:,1), data_filt(:,:,2)]';
    if radius_fun(data_(1:2,end)) > 0.75
        data_(1:2,:) = data_(1:2,:) - repmat(data_(1:2,end), [1 length(data_)]);
        data_(3:4,end) = zeros(2,1);
    end    
    Data = [Data data_];
    clear data_
end

% Position/Velocity Trajectories
Xi_ref     = Data(1:2,:);
Xi_dot_ref = Data(3:end,:);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Step 1: Fit GMM to Trajectory Data        %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GMM Estimation Algorithm %%%%
% 0: Physically-Consistent Non-Parametric (Collapsed Gibbs Sampler)
% 1: GMM-EM Model Selection via BIC
% 2: GMM via Competitive-EM
% 3: CRP-GMM via Collapsed Gibbs Sampler

est_options = [];
est_options.type       = 1;   % GMM Estimation Alorithm Type    
est_options.maxK       = 10;  % Maximum Gaussians for Type 1/2
est_options.do_plots   = 1;   % Plot Estimation Statistics
est_options.adjusts_C  = 0;   % Adjust Sigmas

% Discover Local Models
sample = 3;
[Priors0, Mu0, Sigma0] = discover_local_models(Xi_ref(:,1:sample:end), Xi_dot_ref(:,1:sample:end), est_options);

%% Extract Cluster Labels
est_K      = length(Priors0); 
Priors = Priors0; Mu = Mu0; Sigma = Sigma0;

% Use init from SEDS
% Priors = Priors_0; Mu = Mu_0(1:2,:); Sigma = Sigma_0(1:2,1:2,:);
est_labels =  my_gmm_cluster(Xi_ref, Priors, Mu, Sigma, 'hard', []);

% Visualize Cluster Parameters on Manifold Data
plotGMMParameters( Xi_ref, est_labels, Mu, Sigma);
limits_ = limits + [-0.015 0.015 -0.015 0.015];
axis(limits_)
switch est_options.type   
    case 0
        title('Physically-Consistent Non-Parametric Mixture Model','Interpreter','LaTex', 'FontSize',15);
    case 1        
        title('Best fit GMM with EM-based BIC Model Selection','Interpreter','LaTex', 'FontSize',15);
    case 3
        title('Standard Non-Parametric Mixture Model','Interpreter','LaTex', 'FontSize',15);
end

%%%% Visualize GMM pdf from learnt parameters (for 1D Datasets)
ml_plot_gmm_pdf(Xi_ref, Priors, Mu, Sigma)
clear ds_gmm; ds_gmm.Mu = Mu; ds_gmm.Sigma = Sigma; 
ds_gmm.Priors = Priors; 

% Adjust Covariance Matrices
if est_options.adjusts_C  == 1
    tot_scale_fact = 1; rel_scale_fact = 0.15;
    Sigma = adjust_Covariances(Sigma0, tot_scale_fact, rel_scale_fact);
    ds_gmm.Sigma;
    % Visualize Cluster Parameters on Manifold Data
    plotGMMParameters( Xi_ref, est_labels, Mu, Sigma);
    ml_plot_gmm_pdf(Xi_ref, Priors, Mu, Sigma)  
    limits_ = limits + [-0.015 0.015 -0.015 0.015];
    axis(limits_)
end    
   
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%  Step 2: ESTIMATE CANDIDATE LYAPUNOV FUNCTION PARAMETERS  %%%%%%%%%
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Vxf] = learn_wsaqf(Data, att_g);
P_g = Vxf.P(:,:,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%  Step 3: ESTIMATE SYSTEM DYNAMICS MATRICES  %%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% DS PARAMETER INITIALIZATION %%%%%%%%%%%%%%%
% Type of constraints to estimate Axi+b or LPV
constr_type = 0;      % 0:'convex':     A' + A < 0
                      % 1:'non-convex': A'P + PA < 0
                      % 2:'non-convex': A'P + PA < -Q given P              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%  LPV system sum_{k=1}^{K}h_k(xi)(A_kxi + b_k) %%%%%%%%            
[A_g, b_g, P_est] = optimize_lpv_ds_from_data(Data, att_g, constr_type, ds_gmm, P_g);
if constr_type == 1
    P_g = P_est;
end
        
%%%% FOR DEBUGGING: Check Negative-Definite Constraint %%%%
if constr_type == 0
    constr_violations = zeros(1,est_K);
    for k=1:est_K
        A_t = A_g(:,:,k) + A_g(:,:,k)';
        constr_violations(1,k) = sum(eig(A_t) > 0); % sufficient
    end
    % Check Constraint Violation
    if sum(constr_violations) > 0
        warning(sprintf('Strict System Matrix Constraints are NOT met..'))
    else
        fprintf('All Sufficient System Matrix Constraints are met..\n')
    end
else
    suff_constr_violations = zeros(1,est_K);
    for k=1:est_K
        Pg_A =  A_g(:,:,k)'*P_g + P_g*A_g(:,:,k);
        suff_constr_violations(1,k) = sum(eig(Pg_A + Pg_A') > 0); % strict
    end
    % Check Constraint Violation
    if sum(suff_constr_violations) > 0
        warning(sprintf('Sufficient System Matrix Constraints are NOT met..'))
    else
        fprintf('All Sufficient System Matrix Constraints are met..\n')
    end
end
%% %%%%%%%%%%%%    Plot Resulting DS  %%%%%%%%%%%%%%%%%%%
% Create DS function handle
ds_lpv = @(x) lpv_ds(x, ds_gmm, A_g, b_g);

if ~exist('fig1','var');     fig1 = figure('Color',[1 1 1]);    end
if exist('hs','var');     delete(hs);    end
if exist('hd','var');     delete(hd);    end
[hd] = scatter(Xi_ref(1,:),Xi_ref(2,:),10,[1 0 0],'filled'); hold on
[hs] = plot_ds_model(fig1, ds_lpv, att_g, limits,'medium'); hold on;
limits_ = limits + [-0.015 0.015 -0.015 0.015];
axis(limits_)
box on
grid on
title('GMM LPV-DS $\dot{\xi}=\sum_{k=1}^K\gamma^k(\xi)(A_k(\xi) + b_k)$', 'Interpreter','LaTex','FontSize',20)
xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);

%% Simulate Passive DS Controller function
if with_robot
    struct_stiff = [];
    struct_stiff.DS_type = 'global'; % Options: SEDS, Global-LPV, LAGS, LMDS ?
    struct_stiff.gmm = ds_gmm;
    struct_stiff.A_g = A_g;
    struct_stiff.basis = 'D'; % Options: D (using the same basis as D) or I (using the world basis)
    
    if with_robot
        dt = 0.01;
        simulate_passiveDS(fig1, robot, base, ds_lags_multi, att_g, dt,struct_stiff);
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot WSAQF Lyapunov Function and derivative -- NEW
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Type of plot
contour = 0; % 0: surf, 1: contour
clear lyap_fun_comb lyap_der 

% Lyapunov function
lyap_fun_comb = @(x)lyapunov_function_combined(x, att_g, [], 1, P_g, zeros(2,2), ds_gmm);
title_string = {'$V(\xi) = (\xi-\xi_g^*)^TP_g(\xi-\xi_g^*)$'};

% Derivative of Lyapunov function (gradV*f(x))
lyap_der = @(x)lyapunov_combined_derivative(x, att_g, [], ds_lpv, 1, P_g, zeros(2,2));
title_string_der = {'Lyapunov Function Derivative $\dot{V}_{g}(\xi)$'};

if exist('h_lyap','var');     delete(h_lyap);     end
if exist('h_lyap_der','var'); delete(h_lyap_der); end
h_lyap     = plot_lyap_fct(lyap_fun_comb, contour, limits,  title_string, 0);
h_lyap_der = plot_lyap_fct(lyap_der, contour, limits_,  title_string_der, 1);

%% Compare Velocities from Demonstration vs DS
% Simulated velocities of DS converging to target from starting point
xd_dot = []; xd = [];
% Simulate velocities from same reference trajectory
for i=1:length(Data)
    xd_dot_ = ds_lags_multi(Data(1:2,i));    
    % Record Trajectories
    xd_dot = [xd_dot xd_dot_];        
end

% Plot Demonstrated Velocities vs Generated Velocities
if exist('h_vel','var');     delete(h_vel);    end
h_vel = figure('Color',[1 1 1]);
plot(Data(3,:)', '.-','Color',[0 0 1], 'LineWidth',2); hold on;
plot(Data(4,:)', '.-','Color',[1 0 0], 'LineWidth',2); hold on;
plot(xd_dot(1,:)','--','Color',[0 0 1], 'LineWidth', 1); hold on;
plot(xd_dot(2,:)','--','Color',[1 0 0], 'LineWidth', 1); hold on;
grid on;
legend({'$\dot{\xi}^{ref}_{x}$','$\dot{\xi}^{ref}_{y}$','$\dot{\xi}^{d}_{x}$','$\dot{\xi}^{d}_{y}$'}, 'Interpreter', 'LaTex', 'FontSize', 15)
