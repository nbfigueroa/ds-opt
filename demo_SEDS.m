%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo Code for SEDS Learning Comparison for paper:                       %
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
% Select initialization type
init_type = 'seds-init'; 
% 'seds-init': follows the initialization given in the SEDS code
do_ms_bic = 1;
% 0: Manually set the # of Gaussians
% 1: Do Model Selection with BIC

% 'phys-gmm': provides different initializations
phys_gmm_type = 0;
% 0: Physically-Consistent Non-Parametric (Collapsed Gibbs Sampler)
% 1: GMM-EM Model Selection via BIC
% 2: GMM via Competitive-EM
% 3: CRP-GMM via Collapsed Gibbs Sampler

switch init_type
    case 'seds-init'
        if do_ms_bic
            est_options = [];
            est_options.type       = 1;   % GMM Estimation Alorithm Type
            est_options.maxK       = 10;  % Maximum Gaussians for Type 1/2
            est_options.do_plots   = 1;   % Plot Estimation Statistics
            
            % Discover Local Models
            sample = 1;
            [Priors0, Mu0, Sigma0] = discover_local_models(Xi_ref(:,1:sample:end), Xi_dot_ref(:,1:sample:end), est_options);
            nb_gaussians = length(Priors0);
        else
            % Select manually the number of Gaussian components
            nb_gaussians = 4;
        end
        [Priors_0, Mu_0, Sigma_0] = initialize_SEDS([Xi_ref(:,1:sample:end); Xi_dot_ref(:,1:sample:end)],nb_gaussians); %finding an initial guess for GMM's parameter
        title_string = '$\theta_{\gamma}=\{\pi_k,\mu^k,\Sigma^k\}$ Initial Estimate';    
        
    case 'phys-gmm'
        
        est_options = [];
        est_options.type       = 0;   % GMM Estimation Alorithm Type
        est_options.maxK       = 10;  % Maximum Gaussians for Type 1/2
        est_options.do_plots   = 1;   % Plot Estimation Statistics
        
        % Discover Local Models
        sample = 3;
        [Priors0, Mu0, Sigma0] = discover_local_models(Xi_ref(:,1:sample:end), Xi_dot_ref(:,1:sample:end), est_options);
        nb_gaussians = length(Priors0);
        
        % Using phys-gmm as initialization
        % --- I assume I ran the sampler already in another script        
        Priors_0 = Priors0;
        % Find the corresponding means in the joint-space
        Idx = knnsearch(Xi_ref',Mu0','k',5);
        Mu_0 = zeros(4,nb_gaussians);
        Mu_0(1:2,:) = Mu0;
        for k=1:nb_gaussians
            Mu_0(3:4,k) = mean(Xi_dot_ref(:,Idx(k,:)),2);
        end
        
        d_i =  my_distX2Mu([Xi_ref;Xi_dot_ref], Mu_0, 'L2');
        [~, k_i] = min(d_i, [], 1);
        Sigma_0 = zeros(4,4,nb_gaussians);
        for k=1:nb_gaussians
            idtmp = find(k_i == k);
            % For DS input
            Sigma_0(1:2,1:2,k) = Sigma0(:,:,k);
            
            % For DS output
            Sigma_0(3:4,3:4,k) = cov([Xi_dot_ref(:,idtmp) Xi_dot_ref(:,idtmp)]');            
            % Add a tiny variance to avoid numerical instability
            Sigma_0(3:4,3:4,k) = Sigma_0(3:4,3:4,k) + 1E-5.*diag(ones(2,1));
        end
%         [Priors_0, Mu_0(3:4,:), Sigma_0(3:4,3:4,:)] = EM(data_0(3:4,:), Priors_0, Mu_0(3:4,:), Sigma_0(3:4,3:4,:));
        [Priors_0, Mu_0, Sigma_0] = EM([Xi_ref;Xi_dot_ref], Priors_0, Mu_0, Sigma_0);
        title_string = 'GMM Parameters Initialized with Phys-GMM';    
end
 
% Plot Initial Estimate 
figure('Color', [1 1 1]);
est_labels =  my_gmm_cluster(Xi_ref, Priors_0, Mu_0(1:2,:), Sigma_0(1:2,1:2,:), 'hard', []);
plotGMMParameters( Xi_ref, est_labels, Mu_0(1:2,:), Sigma_0(1:2,1:2,:),1);
title(title_string, 'Interpreter', 'LaTex', 'FontSize',20)
limits_ = limits + [-0.015 0.015 -0.015 0.015];
axis(limits_)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Step 2: Run SEDS Solver        %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear options;
options.tol_mat_bias = 10^-6; % A very small positive scalar to avoid
                              % instabilities in Gaussian kernel [default: 10^-1]                             
options.display = 1;          % An option to control whether the algorithm
                              % displays the output of each iterations [default: true]                            
options.tol_stopping=10^-9;   % A small positive scalar defining the stoppping
                              % tolerance for the optimization solver [default: 10^-10]
options.max_iter = 500;      % Maximum number of iteration for the solver [default: i_max=1000]


options.objective = 'mse';    % 'likelihood'
[Priors Mu Sigma]= SEDS_Solver(Priors_0,Mu_0,Sigma_0,[Xi_ref(:,1:sample:end); Xi_dot_ref(:,1:sample:end)],options); %running SEDS optimization solver
ds_seds = @(x) GMR_SEDS(Priors,Mu,Sigma,x,1:2,3:4);

% Plot SEDS model
fig3 = figure('Color',[1 1 1]);
scatter(att_g(1,1),att_g(2,1),50,[0 0 0],'filled'); hold on
scatter(Xi_ref(1,:),Xi_ref(2,:),10,[1 0 0],'filled'); hold on
plot_ds_model(fig3, ds_seds, att_g, limits,'medium'); hold on;
limits_ = limits + [-0.015 0.015 -0.015 0.015];
axis(limits_)
box on
grid on
xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);
switch options.objective
    case 'mse'        
        title('SEDS Dynamics with $J(\theta_{\gamma})$=MSE', 'Interpreter','LaTex','FontSize',20)
    case 'likelihood'
        title('SEDS Dynamics with $J(\theta_{\gamma})$= log-Likelihood', 'Interpreter','LaTex','FontSize',20)
end


figure('Color', [1 1 1]);
est_labels =  my_gmm_cluster(Xi_ref, Priors', Mu(1:2,:), Sigma(1:2,1:2,:), 'hard', []);
plotGMMParameters( Xi_ref, est_labels, Mu(1:2,:), Sigma(1:2,1:2,:),1);
scatter(Xi_ref(1,:),Xi_ref(2,:),10,[1 0 0],'filled'); hold on
limits_ = limits + [-0.015 0.015 -0.015 0.015];
axis(limits_)
box on
grid on
title('$\theta_{\gamma}=\{\pi_k,\mu^k,\Sigma^k\}$ after SEDS Optimization', 'Interpreter', 'LaTex','FontSize',20)

%% Compare Velocities from Demonstration vs DS
% Simulated velocities of DS converging to target from starting point
xd_dot = []; xd = [];
% Simulate velocities from same reference trajectory
for i=1:length(Data)
    xd_dot_ = ds_seds(Data(1:2,i));    
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

