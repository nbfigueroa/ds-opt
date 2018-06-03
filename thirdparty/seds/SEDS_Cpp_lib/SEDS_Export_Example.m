%% Reading demonstrations from text files
figure
hold on;grid on
clear demos t
num_demo = 0;
for f = 1:9;
    num_demo = num_demo+1;
    data = dlmread(['demonstrations/recording_' num2str(f) '.txt']);

    x = data(:,9:11); %smoothing data
    plot3(x(:,1),x(:,2),x(:,3),'r.')
    demos{num_demo} = x';
    t{num_demo} = data(:,1)';
end
xlabel('x')
ylabel('y')
zlabel('z')
view(-40,40)
axis equal
title('The original demonstrations from the WAM robot.')

%% Preprocessing Data
tol_cutting = 0.005; % A threshold on velocity that will be used for trimming demos
[x0 , xT, Data, index] = preprocess_demos(demos,t,tol_cutting); %preprocessing data
figure
plot3(Data(1,:),Data(2,:),Data(3,:),'r.')
d = size(Data,1)/2;
xlabel('x')
ylabel('y')
zlabel('z')
view(-44,28)
title('Demonstrations after the preprocessing step.')
axis equal
pause(0.1)
%% SEDS Parameters

K = 4; %Number of Gaussian funcitons

% A set of options that will be passed to the solver. Please type 
% 'doc preprocess_demos' in the MATLAB command window to get detailed
% information about other possible options.
options.tol_mat_bias = 10^-10; % A very small positive scalar to avoid
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

%% SEDS Learning Algorithm
[Priors_0, Mu_0, Sigma_0] = initialize_SEDS(Data,K); %finding an initial guess for GMM's parameter
[Priors Mu Sigma]=SEDS_Solver(Priors_0,Mu_0,Sigma_0,Data,options); %running SEDS optimization solver

%% Simulation

% A set of options that will be passed to the Simulator. Please type 
% 'doc preprocess_demos' in the MATLAB command window to get detailed
% information about each option.
opt_sim.dt = 0.01;
opt_sim.i_max = 3000;
opt_sim.tol = 0.05;
d = size(Data,1)/2; %dimension of data
x0_all = Data(1:d,index(1:end-1)); %finding initial points of all demonstrations
fn_handle = @(x) GMR(Priors,Mu,Sigma,x,1:d,d+1:2*d);
[x xd]=Simulation(x0_all,[],fn_handle,opt_sim); %running the simulator
plot3(Data(1,:),Data(2,:),Data(3,:),'r.')

title('Reproductions from the trained model.')
xlabel('x')
ylabel('y')
zlabel('z')
view(-44,28)
axis equal
%% Exporting the trained model to the SEDS lib
out = export2SEDS_Cpp_lib('mySEDSModel.txt',Priors, Mu, Sigma);