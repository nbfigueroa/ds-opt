function [Priors, Mu, Sigma] = discover_local_models(Xi_ref, Xi_dot_ref, est_options)

% Parse Options
est_type       = est_options.type;
max_gaussians  = est_options.maxK;
do_plots       = est_options.do_plots;
exp_scaling    = est_options.exp_scaling;

if isempty(est_options.fixed_K)
    fixed_K        = 0;
else
    fixed_K = est_options.fixed_K;
end
[M N] = size(Xi_ref);

switch est_type
    case 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Option1: Non-parametric Clustering with Pos-Vel-cos-sim prior %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        
        % Compute estimate of length-scale
        [D, max_hist_D, max_D] = computePairwiseDistances(Xi_ref',1);
        sigma = sqrt(max_hist_D/2)
        l = 1/(2*sigma^2)
%         close all;
        
        % Compute element-wise cosine similarities
        S = zeros(length(Xi_dot_ref),length(Xi_dot_ref));
        for i=1:length(Xi_dot_ref)
            for j=1:length(Xi_dot_ref)
                % Compute Velocity component
                % Normalize vectors
                xi_dot_i = Xi_dot_ref(:,i)/norm(Xi_dot_ref(:,i));
                xi_dot_j = Xi_dot_ref(:,j)/norm(Xi_dot_ref(:,j));
                % Compute Angle
                s_angle = atan2(xi_dot_i(2),xi_dot_i(1))-atan2(xi_dot_j(2),xi_dot_j(1));
                 
                % Compute shifted cosine of angle
                cos_angle = cos(s_angle);
                if isnan(cos_angle)
                    cos_angle = 0;
                end
                s = 1 + cos_angle;
                
                % Compute Position component
                xi_i = Xi_ref(:,i);
                xi_j = Xi_ref(:,j);

                % LASA DATASET
                if exp_scaling                     
                    p = exp(-l*norm(xi_i - xi_j));
                else
                    p = 1;
                end
                
                % With Euclidean pairwise kernel
                f = p * s;
                
                % Exponential decay function
%                 f = 2 * (1-exp(-s));
                
                % Shifted Cosine Similarity of velocity vectors
                S(i,j) = f;
                
            end
        end
        
        % Plot Similarity matrix
        if do_plots
            if exist('h_sim','var');     delete(h_sim);    end;
            title_str = 'Physically-Consistent Similarity Confusion Matrix';
            h_sim = plotSimilarityConfMatrix(S, title_str);
            pause(0);
        end
        
        % Setting sampler/model options (i.e. hyper-parameters, alpha, Covariance matrix)
        Xi_ref_mean = mean(Xi_ref,2);
        options                 = [];
        options.type            = 'full';  % Type of Covariance Matrix: 'full' = NIW or 'Diag' = NIG
        options.T               = 20;     % Sampler Iterations
        options.alpha           = max(0.1,0.1*(randi(11)-2)); % Concentration parameter
        
        % Standard Base Distribution Hyper-parameter setting
        if strcmp(options.type,'diag')
            lambda.alpha_0       = M;                    % G(sigma_k^-1|alpha_0,beta_0): (degrees of freedom)
            lambda.beta_0        = sum(diag(cov(Xi_ref')))/M; % G(sigma_k^-1|alpha_0,beta_0): (precision)
        end
        if strcmp(options.type,'full')
            lambda.nu_0        = M;                           % IW(Sigma_k|Lambda_0,nu_0): (degrees of freedom)
            lambda.Lambda_0    = eye(M)*sum(diag(cov(Xi_ref')))/M; % IW(Sigma_k|Lambda_0,nu_0): (Scale matrix)
            %     lambda.Lambda_0    = diag(diag(cov(Xi_ref')));       % IW(Sigma_k|Lambda_0,nu_0): (Scale matrix)
        end
        lambda.mu_0             = Xi_ref_mean;    % hyper for N(mu_k|mu_0,kappa_0)
        lambda.kappa_0          = 1;            % hyper for N(mu_k|mu_0,kappa_0)       
        
        % Run Collapsed Gibbs Sampler
        options.lambda    = lambda;
        options.verbose   = 1;
        [Psi Psi_Stats]   = run_ddCRP_sampler(Xi_ref, S, options);
        est_labels        = Psi.Z_C';
        
        % Extract Learnt cluster parameters
        unique_labels = unique(est_labels);                                
        est_K = length(unique_labels);
        Priors = zeros(1, est_K);
        for k=1:est_K
            Priors(k) = sum(est_labels==unique_labels(k))/length(Xi_ref);
        end
        Mu     = Psi.Theta.Mu;
        Sigma  = Psi.Theta.Sigma;        
        if do_plots
            if exist('h1b','var') && isvalid(h1b), delete(h1b);end
            stats_options = [];
            stats_options.dataset      = 'Drawn Trajectory Data';
            stats_options.true_labels  = [];
            stats_options.Psi          = Psi;
            [ h1b ] = plotSamplerStats( Psi_Stats, stats_options );
        end
        
    case 1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Option2: Cluster Trajectories with GMM-EM + BIC Model Selection %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if fixed_K == 0
            em_type = 'matlab'; repetitions = 10;
            [bic_scores, k] = fit_gmm_bic([Xi_ref],max_gaussians, repetitions, em_type, do_plots);
        else
            k = fixed_K;
        end
        % Train GMM with Optimal k
        warning('off', 'all'); % there are a lot of really annoying warnings when fitting GMMs
        GMM_full = fitgmdist([Xi_ref]', k, 'Start', 'plus', 'CovarianceType','full', 'Regularize', .000001, 'Replicates', 10); %fit a GMM to our data
        warning('on', 'all');
        
        % Extract Model Parameters
        Priors = GMM_full.ComponentProportion;
        Mu = transpose(GMM_full.mu);
        Sigma = GMM_full.Sigma;                
        
    case 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Option3: Cluster Trajectories with GMM using Modified Competitive EM Algorithm (MixEst) %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimation options
        options.verbosity = 1;
        options.sm.numMax = 10;
        options.inner.partial.maxiter = 5;
        options.solver = 'lbfgs'; % options: 'lbfgs' or 'cg'
        options.tolCostDiff = 1e-3;
        options.maxiter = 50;
        
        % Visualization and plotting options
        if do_plots
            figure('Units', 'normalized', 'OuterPosition', [0 0 0.35 0.35], 'Color', [1 1 1])
            options.visualization.axes = subplot(2,2,[1 3]);
            options.plotCost.axes = subplot(2,2,2);
            options.plotGradNorm.axes = subplot(2,2,4);
            options.visualization.mixture.colorize = true;
            options.visualization.stopOnClose = true;
        end
        
        % Run Competitive EM Optimization
        [theta, D] = cem(Xi_ref, options);
        
        % Estimated Parameters
        K = length(theta.p);
        Priors = theta.p';
        Mu     = zeros(M, K);
        Sigma  = zeros(M, M, K);
        for k=1:K
            Mu(:,k)    = theta.D{k}.mu;
            Sigma(:,:, k) = theta.D{k}.sigma;
        end      
        
        
    case 3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Option4: Cluster Trajectories with Chinese Restaurant Process MM sampler (CRP-GMM) %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % CRP-GMM (Frank-Wood's implementation) -- faster
        iterations = 200;
        [class_id, mean_record, covariance_record, K_record, lP_record, alpha_record] = sampler(Xi_ref, iterations);
        [val , Maxiter]  = max(lP_record);
        est_labels       = class_id(:,Maxiter);
        
        % Visualization and plotting options
        if do_plots
            figure('Color',[1 1 1])
            subplot(2,1,1)
            semilogx(1:iterations, lP_record'); hold on;
            semilogx(Maxiter,lP_record(Maxiter),'ko','MarkerSize',10);
            grid on
            xlabel('Gibbs Iteration','Interpreter','LaTex','Fontsize',20); ylabel('LogPr','Interpreter','LaTex','Fontsize',20)
            xlim([1 iterations])
            legend({'$p(Z|Y, \alpha, \lambda)$'},'Interpreter','LaTex','Fontsize',14)
            title(sprintf('CRP-GMM Sampling results, optimal K=%d at iter=%d', length(unique(est_labels)), Maxiter), 'Interpreter','LaTex','Fontsize',20)            
            subplot(2,1,2)
            stairs(K_record, 'LineWidth',2);
            set(gca, 'XScale', 'log')
            xlim([1 iterations])
            xlabel('Gibbs Iteration','Interpreter','LaTex','Fontsize',20); ylabel('$\Psi$ = Estimated K','Interpreter','LaTex','Fontsize',20);
        end
        % Extract Learnt cluster parameters
        unique_labels = unique(est_labels);
        est_K = length(unique_labels);
        Priors = zeros(1, est_K);
        for k=1:est_K
            Priors(k)    = sum(est_labels==unique_labels(k))/length(Xi_ref(:,1:end));
            
        end
        Mu    = mean_record {Maxiter};
        Sigma = covariance_record{Maxiter};
        
end


end

