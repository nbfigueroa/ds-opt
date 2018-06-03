%% Testing DS modeling approaches

% Prepare Data for DS Learning
subs = 5;
data =   [X_demos(1:subs:end), Z_demos(1:subs:end), U_demos(1:subs:end), W_demos(1:subs:end)]';
data_0 = [X_demos_0(1:subs:end), Z_demos_0(1:subs:end), U_demos_0(1:subs:end), W_demos_0(1:subs:end)]';
target = data(1:2,end);


%% 1) Fixed Original Linear Dynamics
figure('Color',[1 1 1]);
fig1 = subplot(1,3,1);
scatter(target(1,1),target(2,1),50,[0 0 0],'filled'); hold on
A = -[10 0;0 10]; b = [0 0]';
ds_lin = @(x) lin_ds(A,b,x);
scatter(data(1,:),data(2,:),10,[0 0 0],'filled'); hold on
plot_ds_model(fig1, ds_lin, target); hold on;
axis tight
title('Original Linear Dynamics', 'Interpreter','LaTex')

xLimits = get(gca,'XLim');  %# Get the range of the x axis
yLimits = get(gca,'YLim');  %# Get the range of the y axis
limits = [xLimits yLimits];

%% 2) Mix Stable Linear Dynamics approx from Demonstrations
figure('Color',[1 1 1])
n_comp = 9;
em_iterations = 2;
clear options;
options.n_iter = em_iterations;        % Max number of EM iterations
options.solver = 'sedumi';             % Solver
options.criterion = 'mse';             % Solver
options.c_reg = 1e-6;                  % Pos def eps margin
options.c_reg_inv = 5e-1;
options.verbose = 1;                   % Verbose (0-5)
options.warning = true;                % Display warning information

% Prior for the attractor
options.prior.mu = target;
options.prior.sigma_inv = [1 0; 0 1];
lambda = em_mix_lds_inv_max(data, n_comp, options);

% Plot result
scatter(target(1,1),target(2,1),50,[0 0 0],'filled'); hold on
scatter(data(1,:),data(2,:),10,[1 0 0],'filled'); hold on
plot_streamlines_mix_lds(lambda, limits);
axis tight
title('Stable Mix LPV Systems from Demo', 'Interpreter','LaTex')

%% 3) Stable Non-Linear Dynamics approx from Demonstrations
data_0 = Data; data = Data;
target = att_g;
% learn SEDS model
clear options;
options.tol_mat_bias = 10^-6; % A very small positive scalar to avoid
                              % instabilities in Gaussian kernel [default: 10^-1]                             
options.display = 1;          % An option to control whether the algorithm
                              % displays the output of each iterations [default: true]                            
options.tol_stopping=10^-9;   % A small positive scalar defining the stoppping
                              % tolerance for the optimization solver [default: 10^-10]
options.max_iter = 500;      % Maximum number of iteration for the solver [default: i_max=1000]

init_type = 'seds-init';

switch init_type
    case 'seds-init'
        
        est_options = [];
        est_options.type       = 1;   % GMM Estimation Alorithm Type
        est_options.maxK       = 10;  % Maximum Gaussians for Type 1/2
        est_options.do_plots   = 1;   % Plot Estimation Statistics
        
        % Discover Local Models
        sample = 2;
        [Priors0, Mu0, Sigma0] = discover_local_models(Xi_ref(:,1:sample:end), Xi_dot_ref(:,1:sample:end), est_options);
        nb_gaussians = length(Priors0);
        
        % Standard K-Means based Initialization
        [Priors_0, Mu_0, Sigma_0] = initialize_SEDS(data_0,nb_gaussians); %finding an initial guess for GMM's parameter
        title_string = 'GMM Parameters Initialized with Standard EM';    
        
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
        Idx = knnsearch(data_0(1:2,:)',Mu0','k',5);
        Mu0
        Mu_0 = zeros(4,nb_gaussians);
        Mu_0(1:2,:) = Mu0;
        for k=1:nb_gaussians
            Mu_0(3:4,k) = mean(data(3:end,Idx(k,:)),2);
        end
        
        d_i =  my_distX2Mu(data_0, Mu_0, 'L2');
        [~, k_i] = min(d_i, [], 1);
        Sigma_0 = zeros(4,4,nb_gaussians);
        for k=1:nb_gaussians
            idtmp = find(k_i == k);
            % For DS input
            Sigma_0(1:2,1:2,k) = Sigma0(:,:,k);
            
            % For DS output
            Sigma_0(3:4,3:4,k) = cov([data_0(3:4,idtmp) data_0(3:4,idtmp)]');            
            % Add a tiny variance to avoid numerical instability
            Sigma_0(3:4,3:4,k) = Sigma_0(3:4,3:4,k) + 1E-5.*diag(ones(2,1));
        end
        [Priors_0, Mu_0(3:4,:), Sigma_0(3:4,3:4,:)] = EM(data_0(3:4,:), Priors_0, Mu_0(3:4,:), Sigma_0(3:4,3:4,:));
%         [Priors_0, Mu_0, Sigma_0] = EM(data_0, Priors_0, Mu_0, Sigma_0);
        title_string = 'GMM Parameters Initialized with Phys-GMM';    
end
       
figure('Color', [1 1 1]);
est_labels =  my_gmm_cluster(Xi_ref, Priors_0, Mu_0(1:2,:), Sigma_0(1:2,1:2,:), 'hard', []);
plotGMMParameters( Xi_ref, est_labels, Mu_0(1:2,:), Sigma_0(1:2,1:2,:),1);
title(title_string, 'Interpreter', 'LaTex')

%% Run SEDS-Solver
options.objective = 'likelihood';    % 'likelihood'
[Priors Mu Sigma]=SEDS_Solver(Priors_0,Mu_0,Sigma_0,data_0,options); %running SEDS optimization solver
ds_seds = @(x) GMR_SEDS(Priors,Mu,Sigma,x,1:2,3:4);

% Plot SEDS model
fig3 = figure('Color',[1 1 1]);
scatter(target(1,1),target(2,1),50,[0 0 0],'filled'); hold on
scatter(data(1,:),data(2,:),10,[1 0 0],'filled'); hold on
plot_ds_model(fig3, ds_seds, att_g, limits,'medium'); hold on;
axis tight
switch options.objective
    case 'mse'        
        title('SEDS Dynamics with Minizimized MSE', 'Interpreter','LaTex')
    case 'likelihood'
        title('SEDS Dynamics with Maximized Likelihood', 'Interpreter','LaTex')
end
%%
figure('Color', [1 1 1]);
est_labels =  my_gmm_cluster(Xi_ref, Priors', Mu(1:2,:), Sigma(1:2,1:2,:), 'hard', []);
plotGMMParameters( Xi_ref, est_labels, Mu(1:2,:), Sigma(1:2,1:2,:),1);
title('GMM Parameters After SEDS Solver', 'Interpreter', 'LaTex')
