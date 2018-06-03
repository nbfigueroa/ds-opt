% This is a matlab script illustrating how to use SEDS_lib to learn
% an arbitrary model from a set of demonstrations.
%%
% To run this demo you need to provide the variable demos composed of all
% demosntration trajectories. To get more detailed information about the
% structure of the variable 'demo', type 'doc preprocess_demos' in the
% MATLAB command window

%% User Parameters and Setting
load('models/recorded_motions/GShape','demos')
% the variable 'demos' composed of 3 demosntrations. Each demonstrations is
% recorded from Tablet-PC at 50Hz. Datas are in millimeters.

% Pre-processing
dt = 0.1; %The time step of the demonstrations
tol_cutting = 1; % A threshold on velocity that will be used for trimming demos

% Training parameters
K = 6; %Number of Gaussian funcitons

% A set of options that will be passed to the solver. Please type 
% 'doc preprocess_demos' in the MATLAB command window to get detailed
% information about other possible options.
options.tol_mat_bias = 10^-6; % A very small positive scalar to avoid
                              % instabilities in Gaussian kernel [default: 10^-15]
                              
options.display = 1;          % An option to control whether the algorithm
                              % displays the output of each iterations [default: true]
                              
options.tol_stopping=10^-10;  % A small positive scalar defining the stoppping
                              % tolerance for the optimization solver [default: 10^-10]

options.max_iter = 500;       % Maximum number of iteration for the solver [default: i_max=1000]

options.objective = 'mse';    % 'likelihood': use likelihood as criterion to
                              % optimize parameters of GMM
                              % 'mse': use mean square error as criterion to
                              % optimize parameters of GMM
                              % 'direction': minimize the angle between the
                              % estimations and demonstrations (the velocity part)
                              % to optimize parameters of GMM                              
                              % [default: 'mse']

%% Putting GMR and SEDS library in the MATLAB Path
if isempty(regexp(path,['SEDS_lib' pathsep], 'once'))
    addpath([pwd, '/SEDS_lib']);    % add SEDS dir to path
end
if isempty(regexp(path,['GMR_lib' pathsep], 'once'))
    addpath([pwd, '/GMR_lib']);    % add GMR dir to path
end

%% SEDS learning algorithm
[x0 , xT, Data, index] = preprocess_demos(demos,dt,tol_cutting); %preprocessing datas
[Priors_0, Mu_0, Sigma_0] = initialize_SEDS(Data,K); %finding an initial guess for GMM's parameter
[Priors Mu Sigma]=SEDS_Solver(Priors_0,Mu_0,Sigma_0,Data,options); %running SEDS optimization solver

%% Simulation

% A set of options that will be passed to the Simulator. Please type 
% 'doc preprocess_demos' in the MATLAB command window to get detailed
% information about each option.
opt_sim.dt = 0.1;
opt_sim.i_max = 3000;
opt_sim.tol = 0.1;
d = size(Data,1)/2; %dimension of data
x0_all = Data(1:d,index(1:end-1)); %finding initial points of all demonstrations
fn_handle = @(x) GMR(Priors,Mu,Sigma,x,1:d,d+1:2*d);
[x xd]=Simulation(x0_all,[],fn_handle,opt_sim); %running the simulator
%%
% plotting the result
figure('name','Results from Simulation','position',[265   200   520   720])
sp(1)=subplot(3,1,1);
hold on; box on
% plotGMM(Mu(1:2,:), Sigma(1:2,1:2,:), [0.6 1.0 0.6], 1,[0.6 1.0 0.6]);
plot(Data(1,:),Data(2,:),'r.')
xlabel('$\xi_1 (mm)$','interpreter','latex','fontsize',15);
ylabel('$\xi_2 (mm)$','interpreter','latex','fontsize',15);
title('Simulation Results')

sp(2)=subplot(3,1,2);
hold on; box on
% plotGMM(Mu([1 3],:), Sigma([1 3],[1 3],:), [0.6 1.0 0.6], 1,[0.6 1.0 0.6]);
plot(Data(1,:),Data(3,:),'r.')
xlabel('$\xi_1 (mm)$','interpreter','latex','fontsize',15);
ylabel('$\dot{\xi}_1 (mm/s)$','interpreter','latex','fontsize',15);

sp(3)=subplot(3,1,3);
hold on; box on
% plotGMM(Mu([2 4],:), Sigma([2 4],[2 4],:), [0.6 1.0 0.6], 1,[0.6 1.0 0.6]);
plot(Data(2,:),Data(4,:),'r.')
xlabel('$\xi_2 (mm)$','interpreter','latex','fontsize',15);
ylabel('$\dot{\xi}_2 (mm/s)$','interpreter','latex','fontsize',15);

for i=1:size(x,3)
    plot(sp(1),x(1,:,i),x(2,:,i),'linewidth',2)
    plot(sp(2),x(1,:,i),xd(1,:,i),'linewidth',2)
    plot(sp(3),x(2,:,i),xd(2,:,i),'linewidth',2)
    plot(sp(1),x(1,1,i),x(2,1,i),'ok','markersize',5,'linewidth',5)
    plot(sp(2),x(1,1,i),xd(1,1,i),'ok','markersize',5,'linewidth',5)
    plot(sp(3),x(2,1,i),xd(2,1,i),'ok','markersize',5,'linewidth',5)
end

for i=1:3
    axis(sp(i),'tight')
    ax=get(sp(i));
    axis(sp(i),...
        [ax.XLim(1)-(ax.XLim(2)-ax.XLim(1))/10 ax.XLim(2)+(ax.XLim(2)-ax.XLim(1))/10 ...
        ax.YLim(1)-(ax.YLim(2)-ax.YLim(1))/10 ax.YLim(2)+(ax.YLim(2)-ax.YLim(1))/10]);
    plot(sp(i),0,0,'k*','markersize',15,'linewidth',3)
    if i==1
        D = axis(sp(i));
    end
end

% plotting streamlines
figure('name','Streamlines','position',[800   90   560   320])
plotStreamLines(Priors,Mu,Sigma,D)
hold on
plot(Data(1,:),Data(2,:),'r.')
plot(0,0,'k*','markersize',15,'linewidth',3)
xlabel('$\xi_1 (mm)$','interpreter','latex','fontsize',15);
ylabel('$\xi_2 (mm)$','interpreter','latex','fontsize',15);
title('Streamlines of the model')
set(gca,'position',[0.1300    0.1444    0.7750    0.7619])