%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Demo for testing non-linear DS estimation via diff. matching from Perrin %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 1 (DATA LOADING): Load Datasets %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc
%%%%%%%%%%%%%%%%%%%%%%%%% Select a Dataset %%%%%%%%%%%%%%%%%%%%%%%%%%
% 1:     Messy Snake Dataset   (2D)
% 2:     L-shape Dataset       (2D)
% 3:     A-shape Dataset       (2D)
% 4:     S-shape Dataset       (2D)
% (N/A)  Dual-behavior Dataset (2D)
% 6:     Via-point Dataset     (3D) -- 15 trajectories recorded at 100Hz
% 7:     Sink Dataset          (3D) -- 21 trajectories recorded at 100Hz
% 8:     CShape top            (3D) -- 10 trajectories recorded at 100Hz
% 9:     CShape bottom         (3D) -- 10 trajectories recorded at 100Hz
% (N/A)  CShape all            (3D) -- 20 trajectories recorded at 100Hz
% (N/A)  Cube arranging        (3D) -- 20 trajectories recorded at 100Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg_dir         = '/home/nbfigueroa/Dropbox/PhD_papers/CoRL-2018/code/ds-opt/';
chosen_dataset  = 2; 
sub_sample      = 1; % '>2' for real 3D Datasets, '1' for 2D toy datasets
nb_trajectories = 0; % For real 3D data
[Data, Data_sh, att, x0_all, data, dt] = load_dataset_DS(pkg_dir, chosen_dataset, sub_sample, nb_trajectories);

% Position/Velocity Trajectories
vel_samples = 10; vel_size = 0.5; 
[h_data, h_att, h_vel] = plot_reference_trajectories_DS(Data, att, vel_samples, vel_size);

% Some options
step_demos  = 1;  %Step the points
doNormalize = 0;  %Normalize the demos using the variance in each direction; This works better in general
nb_demos    = length(data);

%Options concerning the transformation (NOTE:: Normalization deforms the trajectories!)
if doNormalize
    optsSearch = {'maxCoef', '10', 'nb_iteration','150', 'regularise', '5e-4', 'conv_crit', '1e-6', 'division_coefList', '[3,3,2.5,2.0,1.5,1.1,1.1]', 'safeCoeffList', '0.6*[1, 1, 1, 1, 1, 1, 1]', 'doPlot', '0'};
else
    optsSearch = {'maxCoef', '5', 'nb_iteration','150', 'regularise', '1e-3', 'conv_crit', '1e-6', 'division_coefList', '[3,3,2.5,2.0,1.5,1.1,1.1]', 'safeCoeffList', '0.55*[1, 1, 1, 1, 1, 1, 1]', 'doPlot', '0'};
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2: Estimate Diffeomorphic Matching Function Parameters  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
source = [linspace(target_trajectory(1,1), target_trajectory(1,end), lX); linspace(target_trajectory(2,1), target_trajectory(2,end), lX)];
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

%Max for plotting
XLimPlot = [min([1.2*alltarget_trajectory, 0.8*alltarget_trajectory, -.25*ones(dim,1)],[], 2), max([1.2*alltarget_trajectory, 0.8*alltarget_trajectory, .25*ones(dim,1)],[],2)];
XLimPlot2 = [min([1.2*allSource, 0.8*allSource, -.25*ones(dim,1)],[], 2), max([1.2*allSource, 0.8*allSource, .25*ones(dim,1)],[],2)];

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
fig1 = figure('Color',[1 1 1]); hold all;
% Demonstrated Mean Trajectory \xi
scatter(target_trajectory(1,:),target_trajectory(2,:), 20, [1 0 0],'filled');
if nDemos > 1
    % All Demonstrated Trajectories \xi
    scatter(alltarget_trajectory(1,:),alltarget_trajectory(2,:), 10, [0 0 0],'filled');
end
% Virtual Trajectory \chi
scatter(source(1,:),source(2,:), 20, [0 0 1], 'filled'); hold on;

% Demonstrated Trajectory transformed to Virtual through \phi^-1(\xi)
invPhiTarg = inverse_diff_fun(target_trajectory);%Applying the inverse transformation to the mean target_trajectory
scatter(invPhiTarg(1,:), invPhiTarg(2,:), 20, [0 0.5 1], '+'); hold on;

% Virtual Trajectory transformed to Demonstration through \\phi(\chi)
PhiSrc = diff_fun(source);
scatter(PhiSrc(1,:), PhiSrc(2,:), 20, [1 0.5 0], '+'); hold on;

% Plotting the deformed grid through the diffeomorphism
plotGrid(1, XLimPlot2, 10, 1000, diff_fun);
if nDemos > 1    
legend({'Mean Demo trajectory $\Xi=\{\xi_1,\dots,\xi_T\}$','Demonstrated trajectories $\Xi=\{\xi_1,\dots,\xi_T\}$', ... , 
    'Virtual Linear Trajectory $\chi=\{\chi_1,\dots,\chi_T\}$', ...    
    'Transformed Demonstrated Trajectory $\Phi^{-1}(\Xi)=\{\phi^{-1}(\xi_1),\dots,\phi^{-1}(\xi_T)\}$',...
    'Transformed Virtual Trajectory $\Phi(\chi)=\{\phi(\chi_1),\dots,\phi(\chi_T)\}$'},'Interpreter','LaTex','FontSize',10)
else
    legend({'Mean Demo trajectory $\Xi=\{\xi_1,\dots,\xi_T\}$', ... , 
    'Virtual Linear Trajectory $\chi=\{\chi_1,\dots,\chi_T\}$', ...    
    'Transformed Demonstrated Trajectory $\Phi^{-1}(\Xi)=\{\phi^{-1}(\xi_1),\dots,\phi^{-1}(\xi_T)\}$',...
    'Transformed Virtual Trajectory $\Phi(\chi)=\{\phi(\chi_1),\dots,\phi(\chi_T)\}$'},'Interpreter','LaTex','FontSize',10)
end
grid on;
title('Results of Diffeomorphic Matching Algorithm','Interpreter','LaTex', 'FontSize',15)
xlim(XLimPlot(1,:));
ylim(XLimPlot(2,:));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Step 3: Generated Deformed Dynamics function       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate DS function
EIG0 = -diag([1,2.*ones(1,dim-1)]);
ds_diff = @(x) diffeomorphic_ds(x-repmat(att,[1 size(x,2)]), EIG0, source, jac_inverse_diff_fun);

% Plot DS and demonstrations
plot_repr = 1;

fig2 = figure('Color',[1 1 1]);
[hd] = scatter(Data(1,:),Data(2,:),10,[1 0 0],'filled'); hold on
limits = axis;
limits_ = limits + [-0.15 0.15 -0.15 0.15];
[hs] = plot_ds_model(fig1, ds_diff, [0 0]', limits_,'medium'); hold on;
axis(limits_)
box on
grid on
xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);

% Simulate trajectories and plot them on top
if plot_repr
    opt_sim = [];
    opt_sim.dt = 0.01;
    opt_sim.i_max = 3000;
    opt_sim.tol = 0.1;
    opt_sim.plot = 0;
    [x_sim xd_sim]=Simulation(x0_all ,[],ds_diff, opt_sim);
    [hr] = scatter(x_sim(1,:),x_sim(2,:),10,[0 0 0],'filled'); hold on
end
title('Diffeomorphic Dynamics $\dot{\xi} = A(\phi^{-1}(\xi))J_{\phi}(\phi^{-1}(\xi))\phi^{-1}(\xi) $', 'Interpreter','LaTex','FontSize',15)

%% Compare Velocities from Demonstration vs DS
% Simulated velocities of DS converging to target from starting point
xd_dot = []; xd = [];
% Simulate velocities from same reference trajectory
for i=1:length(Data_sh)
    x_ = Data_sh(1:2,i);
    xd_dot_ = feval(ds_diff, x_);
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
legend({'$\dot{\xi}^{ref}_{1}$','$\dot{\xi}^{ref}_{2}$','$\dot{\xi}^{d}_{1}$','$\dot{\xi}^{d}_{2}$'}, 'Interpreter', 'LaTex', 'FontSize', 15)
%% Compute Errors
% Compute RMSE on training data
rmse = mean(rmse_error(ds_lpv, Xi_ref, Xi_dot_ref));
fprintf('LPV-DS with (O%d), got prediction RMSE on training set: %d \n', constr_type+1, rmse);

% Compute e_dot on training data
edot = mean(edot_error(ds_lpv, Xi_ref, Xi_dot_ref));
fprintf('LPV-DS with (O%d), got e_dot on training set: %d \n', constr_type+1, edot);

% Compute DTWD between train trajectories and reproductions
nb_traj       = size(x_sim,3);
ref_traj_leng = size(Xi_ref,2)/nb_traj;
dtwd = zeros(1,nb_traj);
for n=1:nb_traj
    start_id = 1+(n-1)*ref_traj_leng;
    end_id   = n*ref_traj_leng;
   dtwd(1,n) = dtw(x_sim(:,:,n)',Xi_ref(:,start_id:end_id)',20);
end
fprintf('LPV-DS got DTWD of reproduced trajectories: %2.4f +/- %2.4f \n', mean(dtwd),std(dtwd));