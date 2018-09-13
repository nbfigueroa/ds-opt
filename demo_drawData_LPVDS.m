%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo Code for GMM-based LPV_DS Learning for paper:                      %
%  'A Physically-Consistent Bayesian Non-Parametric Mixture Model for     %
%   Dynamical System Learning.'                                           %
% With this script you can load 2D toy trajectories or even real-world 
% trajectories acquired via kinesthetic taching and test the different    %
% GMM fitting approaches.                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA LOADING OPTION 1:  Draw with GUI  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc

fig1 = figure('Color',[1 1 1]);
% Axis limits
limits = [-6 0.5 -0.5 2];
axis(limits)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.55, 0.2646 0.4358]);
grid on

% Draw Reference Trajectories
[data, hp] = draw_mouse_data_on_DS(fig1, limits);
Data = []; x0_all = []; x0_end = []; Data_sh = [];
for l=1:length(data)    
    data_ = data{l};
    x0_end = [x0_end data_(1:2,end)];
    Data = [Data data_];
    x0_all = [x0_all data_(1:2,1)];
    
    % Shift data to origin for (O2)
    data_(1:2,:) = data_(1:2,:) - repmat(data_(1:2,end), [1 length(data_)]);
    data_(3:4,end) = zeros(2,1);

    Data_sh = [Data_sh data_];
end

% Global Attractor of DS
att_g = mean(x0_end,2);

%% Position/Velocity Trajectories
delete(hp)
scatter(att_g(1),att_g(2),100,[0 0 0],'d'); hold on;
if exist('h_att_g','var');  delete(h_att_g); end
[h_data, h_att, h_vel] = plot_reference_trajectories(Data, att_g, [], 10);
grid on;
box on;
title('Demonstrated Trajectories','Interpreter','LaTex','FontSize',20);
xlabel('$x_1$','Interpreter','LaTex','FontSize',20);
ylabel('$x_2$','Interpreter','LaTex','FontSize',20);

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
% 3: CRP-GMM (Collapsed Gibbs Sampler)

est_options = [];
est_options.type        = 0;   % GMM Estimation Alorithm Type    
est_options.maxK        = 15;  % Maximum Gaussians for Type 1/2
est_options.do_plots    = 1;   % Plot Estimation Statistics
est_options.adjusts_C   = 1;   % Adjust Sigmas
est_options.fixed_K     = [];  % Fix K and estimate with EM
est_options.exp_scaling = 1;   % Scaling for the similarity to improve locality

% Discover Local Models
sample = 4;
[Priors0, Mu0, Sigma0] = discover_local_models(Xi_ref(:,1:sample:end), Xi_dot_ref(:,1:sample:end), est_options);

%% Extract Cluster Labels
est_K      = length(Priors0); 
Priors = Priors0; Mu = Mu0; Sigma = Sigma0;
[~, est_labels] =  my_gmm_cluster(Xi_ref, Priors, Mu, Sigma, 'hard', []);

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
        title('Bayesian Non-Parametric Mixture Model','Interpreter','LaTex', 'FontSize',15);
end

%%% Visualize GMM pdf from learnt parameters (for 2D Datasets)
ml_plot_gmm_pdf(Xi_ref, Priors, Mu, Sigma, limits)
clear ds_gmm; ds_gmm.Mu = Mu; ds_gmm.Sigma = Sigma; 
ds_gmm.Priors = Priors; 

% Adjust Covariance Matrices
est_options.adjusts_C  = 1;
if est_options.adjusts_C  == 1
    tot_scale_fact = 1; rel_scale_fact = 0.15;
    Sigma = adjust_Covariances(Sigma0, tot_scale_fact, rel_scale_fact);
    ds_gmm.Sigma = Sigma;
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
        title('Bayesian Non-Parametric Mixture Model','Interpreter','LaTex', 'FontSize',15);
    end   
    
    ml_plot_gmm_pdf(Xi_ref, Priors, Mu, Sigma, limits)  

end    
   
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%  Step 2: ESTIMATE SYSTEM DYNAMICS MATRICES  %%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% DS OPTIMIZATION OPTIONS %%%%%%%%%%%%%%%%%%%%
% Type of constraints/optimization 
constr_type = 0;      % 0:'convex':     A' + A < 0
                      % 1:'non-convex': A'P + PA < 0
                      % 2:'non-convex': A'P + PA < -Q given P                                  
init_cvx    = 0;      % 0/1: initialize non-cvx problem with cvx                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if constr_type == 0 || constr_type == 1
    P_opt = [];
else
    [Vxf] = learn_wsaqf(Data, att_g);
    P_opt = Vxf.P(:,:,1);
end

%%%%%%%%  LPV system sum_{k=1}^{K}\gamma_k(xi)(A_kxi + b_k) %%%%%%%%            
[A_g, b_g, P_est] = optimize_lpv_ds_from_data(Data, att_g, constr_type, ds_gmm, P_opt, init_cvx);

%% %%%%%%%%%%%%    Plot Resulting DS  %%%%%%%%%%%%%%%%%%%
% Create DS function handle
ds_lpv = @(x) lpv_ds(x, ds_gmm, A_g, b_g);

fig1 = figure('Color',[1 1 1]);
% if ~exist('fig1','var');     fig1 = figure('Color',[1 1 1]);    end
% if exist('hs','var');     delete(hs);    end
% if exist('hd','var');     delete(hd);    end
[hd] = scatter(Xi_ref(1,:),Xi_ref(2,:),10,[1 0 0],'filled'); hold on
[hs] = plot_ds_model(fig1, ds_lpv, att_g, limits,'medium'); hold on;
limits_ = limits + [-0.015 0.015 -0.015 0.015];
axis(limits_)
box on
grid on

switch constr_type
    case 0
        title('GMM LPV-DS with (O1) Optimization', 'Interpreter','LaTex','FontSize',20)
    case 1
        title('GMM LPV-DS with (O2) Optimization', 'Interpreter','LaTex','FontSize',20)
    case 2
        title('GMM LPV-DS with (O3) Optimization', 'Interpreter','LaTex','FontSize',20)
end

xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);

% Simulate trajectories and plot them on top
plot_repr = 1;
if plot_repr
    opt_sim = [];
    opt_sim.dt = 0.01;
    opt_sim.i_max = 3000;
    opt_sim.tol = 0.1;
    opt_sim.plot = 0;
    [x_seds xd_seds]=Simulation(x0_all ,[],ds_lpv, opt_sim);
%     if exist('hr','var');     delete(hr);    end
    [hr] = scatter(x_seds(1,:),x_seds(2,:),10,[0 0 0],'filled'); hold on
end


% Compute RMSE on training data
rmse = mean(rmse_error(ds_lpv, Xi_ref, Xi_dot_ref));
fprintf('LPV-DS with (O%d), got prediction RMSE on training set: %d \n', constr_type+1, rmse);

% Compute e_dot on training data
edot = mean(edot_error(ds_lpv, Xi_ref, Xi_dot_ref));
fprintf('LPV-DS with (O%d), got e_dot on training set: %d \n', constr_type+1, edot);

% Compute DTWD between train trajectories and reproductions
nb_traj       = size(x_seds,3);
ref_traj_leng = size(Xi_ref,2)/nb_traj;
dtwd = zeros(1,nb_traj);
for n=1:nb_traj
    start_id = 1+(n-1)*ref_traj_leng;
    end_id   = n*ref_traj_leng;
   dtwd(1,n) = dtw(x_seds(:,:,n)',Xi_ref(:,start_id:end_id)',20);
end
fprintf('LPV-DS got DTWD of reproduced trajectories: %2.4f +/- %2.4f \n', mean(dtwd),std(dtwd));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Plot Choosen Lyapunov Function and derivative  %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Type of plot
contour = 1; % 0: surf, 1: contour
clear lyap_fun_comb lyap_der 

switch constr_type
    case 0 
        P = eye(2);
        title_string = {'$V(\xi) = (\xi-\xi^*)^T(\xi-\xi^*)$'};
    case 1
        P = P_est;
        title_string = {'$V(\xi) = (\xi-\xi^*)^TP(\xi-\xi^*)$'};
    case 2
        P = P_opt;
        title_string = {'$V(\xi) = (\xi-\xi^*)^TP(\xi-\xi^*)$'};
end

% Lyapunov function
lyap_fun = @(x)lyapunov_function_PQLF(x, att_g, P);

% Derivative of Lyapunov function (gradV*f(x))
lyap_der = @(x)lyapunov_derivative_PQLF(x, att_g, P, ds_lpv);
title_string_der = {'Lyapunov Function Derivative $\dot{V}(\xi)$'};

% if exist('h_lyap','var');     delete(h_lyap);     end
% if exist('h_lyap_der','var'); delete(h_lyap_der); end
h_lyap     = plot_lyap_fct(lyap_fun, contour, limits,  title_string, 0);
xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);

h_lyap_der = plot_lyap_fct(lyap_der, contour, limits_,  title_string_der, 1);
[hd] = scatter(Xi_ref(1,:),Xi_ref(2,:),10,[1 1 0],'filled'); hold on
xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);

%% Compare Velocities from Demonstration vs DS
% Simulated velocities of DS converging to target from starting point
xd_dot = []; xd = [];
% Simulate velocities from same reference trajectory
for i=1:length(Data)
    xd_dot_ = ds_lpv(Data(1:2,i));    
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
