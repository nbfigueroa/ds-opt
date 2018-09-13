%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Demo for testing non-linear compliant DS with diff. matching from Perrin %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw Dataset
clear all; close all;
% Draw Reference Trajectory
fig0 = figure('Color', [1 1 1]);
limits = [-3.5 1 -2.5 1];
xlim(limits(1:2));
ylim(limits(3:4));
scatter(0, 0, 50, [0 0.5 1], '+'); hold on;
set(gca,'Ydir','reverse')
data = draw_mouse_data_on_DS(fig0, limits);
%%
data_12 = data{1}(:,1:2);
dt = abs((data_12(1,1) -  data_12(1,2))/data_12(3,1));
dim = 2;
use_demos = []; %Use only the demos with the appearing here; Use all if empty

% Some options
step_demos = 1;  %Step the points
doNormalize = 1; %Normalize the demos using the variance in each direction; This works better in general

%Options concerning the transformation
if doNormalize
    optsSearch = {'maxCoef', '10', 'nb_iteration','150', 'regularise', '5e-4', 'conv_crit', '1e-6', 'division_coefList', '[3,3,2.5,2.0,1.5,1.1,1.1]', 'safeCoeffList', '0.6*[1, 1, 1, 1, 1, 1, 1]', 'doPlot', '0'};
else
    optsSearch = {'maxCoef', '5', 'nb_iteration','150', 'regularise', '1e-3', 'conv_crit', '1e-6', 'division_coefList', '[3,3,2.5,2.0,1.5,1.1,1.1]', 'safeCoeffList', '0.55*[1, 1, 1, 1, 1, 1, 1]', 'doPlot', '0'};
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: Estimate Diffeomorphic Matching Parameters for Given Demonstrations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare data
fprintf('- Preparing Data for Diffeomorphic Matching...');
close all;
if isempty(use_demos)
    use_demos = 1:length(data);
end
[Xinit, target, allTarget, allTargetV, allSource,  ... ,
    allSourceV, nDemos, indD] = processDrawnDataset(data, use_demos, step_demos,dt);
fprintf('done\n');

% Uses the mean of all demonstrations as target Xi
% (The transformation maps the points given in "source" onto the points in "target"
[dim, lX] = size(target);
target = target./nDemos;
targetV = [diff(target,[],2),zeros(dim,1)]./dt;

%Define the source trajectory: A straight line between the initial and
%final point of the mean target trajectory
source = [linspace(target(1,1), target(1,end), lX); linspace(target(2,1), target(2,end), lX)];
sourceV = [diff(source, [], 2).*(step_demos/dt), zeros(dim,1)];

% Sub-sample source/target trajectories depending on the step-size
source  = source(:,indD);
sourceV = sourceV(:,indD);
target  = target(:,indD);
targetV = targetV(:,indD);

% Normalize Trajectories
if doNormalize
    [Xinit, target, targetV, source, sourceV, allTarget, ... ,
        allTargetV, allSource, allSourceV] = normalizeTrajectories(Xinit, target, ... ,
        targetV, source, sourceV, allTarget, allTargetV, allSource, allSourceV);
end

%Max for plotting
XLimPlot = [min([1.2*allTarget, 0.8*allTarget, -.25*ones(dim,1)],[], 2), max([1.2*allTarget, 0.8*allTarget, .25*ones(dim,1)],[],2)];
XLimPlot2 = [min([1.2*allSource, 0.8*allSource, -.25*ones(dim,1)],[], 2), max([1.2*allSource, 0.8*allSource, .25*ones(dim,1)],[],2)];

%Search for the transformation parameters
fprintf('- Searching for Transformation Parameters...');
tic;
[ centers, targets, coefs, division_coef, nb_iteration ] = iterativeSearch( source, target, optsSearch{:} );
toc;
fprintf('done\n')

% Create all lambda functions needed to evaluate the transformation
% Forward transformation (Phi)
diff_fun             = @(pt) result_function(centers, targets, coefs, division_coef, nb_iteration, pt);
% Backward transformation (Phi^-1)
inverse_diff_fun     = @(pt) result_function_reverse(centers, targets, coefs, division_coef, nb_iteration, pt);
% Backward transformation and calculation of the jacobian in point (J^-1)
jac_inverse_diff_fun = @(pt) result_function_reverse_Jac(centers, targets, coefs, division_coef, nb_iteration, pt);

%%%%%%%%%%%% General plot with the transformed sources and targets %%%%%%%%%%%%
if exist('h_trans','var'); delete(h_trans);end
h_trans = figure('Color',[1 1 1]); hold all;

% Demonstrated Mean Trajectory \xi
scatter(target(1,:),target(2,:), 20, [1 0 0],'filled');
if nDemos > 1
    % All Demonstrated Trajectories \xi
    scatter(allTarget(1,:),allTarget(2,:), 10, [0 0 0],'filled');
end
% Virtual Trajectory \chi
scatter(source(1,:),source(2,:), 20, [0 0 1], 'filled'); hold on;

% Demonstrated Trajectory transformed to Virtual through \phi^-1(\xi)
invPhiTarg = inverse_diff_fun(target);%Applying the inverse transformation to the mean target
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
pause(1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Step 2: Plot Streamlines of Deformed Dynamics            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose DS parameters
if exist('hs','var');     delete(hs);    end
if exist('ha1','var');    delete(ha1);   end
if exist('ha2','var');    delete(ha2);   end
if exist('ht','var');     delete(ht);    end
if exist('h_dec','var');  delete(h_dec); end
if exist('h_data','var'); delete(h_data);end
if exist('fig1','var');   delete(fig1);end

fig1 = figure('Color',[1 1 1]);
xlim(XLimPlot(1,:));
ylim(XLimPlot(2,:));
limits = [xlim ylim];

% Demonstrated Trajectory \xi
scatter(PhiSrc(1,:)-att(1),PhiSrc(2,:)-att(2), 10, [1 0 0],'filled'); hold on;

% Plot Dynamical System
ds_diff = @(x) nonlin_diff_ds(x, jac_inverse_diff_fun, target);
fprintf('Plotting DS...');
[hd] = scatter(Xi_ref(1,:)-att(1),Xi_ref(2,:)-att(2),10,[1 0 0],'filled'); hold on
[hs] = plot_ds_model(fig1, ds_diff, [0 0]', limits,'medium'); hold on;
title('Non-linear Diff Dynamics $\dot{\xi} = A(\phi^{-1}(\xi))J_{\phi}(\phi^{-1}(\xi))\phi^{-1}(\xi) $', 'Interpreter','LaTex','FontSize',15)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       OPTIONAL: Plot Resulting Lyapunov function of transformed space       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the resulting Lyapunov functions (transformation applied to
%x²+y²==alpha²
lyapunov_fun = @(x)(sum(x.*x));
virtual_fun = @(x)(x);

if exist('h_lyap','var'); delete(h_lyap);end
h_lyap  = figure( 'Color', [1 1 1]);
subplot(1,2,1)
xlim(XLimPlot(1,:));
ylim(XLimPlot(2,:));
my_LyapFun2D({gcf(), gca()}, lyapunov_fun, virtual_fun,[]); hold on;
% Virtual Trajectory \chi
scatter(source(1,:),source(2,:), 10, [0 0 1], 'filled'); hold on;
% Attractor
scatter(0, 0, 30, [0 0 0], '+'); hold on;
title('Initial Lyapunov Function $V(\chi)$','Interpreter','LaTex', 'FontSize',15);

subplot(1,2,2)
xlim(XLimPlot(1,:));
ylim(XLimPlot(2,:));
my_LyapFun2D({gcf(), gca()}, lyapunov_fun, diff_fun,[]); hold on;
% Demonstrated Trajectory \xi
scatter(target(1,:),target(2,:), 10, [1 0 0],'filled'); hold on;
if nDemos > 1
    % All Demonstrated Trajectories \xi
    scatter(allTarget(1,:),allTarget(2,:), 10, [0 0 0],'filled');
end
% Attractor
scatter(0, 0, 30, [0 0 0], '+'); hold on;

title('Resulting Lyapunov Function $V(\Phi(\xi))$','Interpreter','LaTex', 'FontSize',15);
pause(1);
