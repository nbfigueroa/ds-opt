function demo_Plot_Results
% This is a matlab function illustrating the obtained results using SEDS
% for the 24 handwriting motions provided in the folder 'models'.

%%
% names = {'Angle','Bump','CShape','GShape',...
%          'JShape','JShape_2','Khamesh','Line',...
%          'NShape','PShape','RShape','Saeghe',...
%          'Sharpc','Sine','Soft_Sine','Spoon',...
%          'Sshape','Trapezoid','WShape','Zshape',...
%          'Multi_Models_1','Multi_Models_2','Multi_Models_3','Multi_Models_4'};

names = {'CShape','GShape','JShape_2','Multi_Models_3'};

fprintf('\n\n\nAvailable Models:\n')
for i=1
    for j=1:4
        fprintf('%2u) %-18s',(i-1)*4+j,names{(i-1)*4+j})
    end
    fprintf('\n')
end

n = input('\nType the number of the model you wish to see the result and press enter: ');
if n<1 || n>4
    disp('Wrong model number!')
    disp('Please try again and type a number between 1-4.')
    return
end
%% preprocessing

load(['models/SEDS_models/' names{n}],'Priors','Mu','Sigma','demos','dt','tol_cutting') %loading the model

%<comment>
fprintf('Model is loaded successfully.\n')
%</comment>

if isempty(regexp(path,['SEDS_lib' pathsep], 'once'))
    addpath([pwd, '/SEDS_lib']);    % add SEDS dir to path
end
if isempty(regexp(path,['GMR_lib' pathsep], 'once'))
    addpath([pwd, '/GMR_lib']);    % add GMR dir to path
end

[tmp , tmp, Data, index] = preprocess_demos(demos,dt,tol_cutting); %preprocessing datas

%% Simulation

% A set of options that will be passed to the Simulator. Please type 
% 'doc preprocess_demos' in the MATLAB command window to get detailed
% information about each option.
opt_sim.dt = 0.05;
opt_sim.i_max = 3000;
opt_sim.tol = 0.1;

d = size(Data,1)/2; %dimension of data
x0_all = Data(1:d,index(1:end-1)); %finding initial points of all demonstrations

%<comment>
fprintf('\nNow we can run simulator starting from initial points of all demos.\n')
fprintf('Press Enter to start simulator.\n')
pause
%</comment>
fn_handle = @(x) GMR(Priors,Mu,Sigma,x,1:d,d+1:2*d);
[x, xd, tmp, xT]=Simulation(x0_all,[],fn_handle,opt_sim); %running the simulator
D = plot_results(x,xd,Data,xT,Mu,Sigma);

%<comment>
fprintf('We continue our evaluation by plotting streamlines of the model.\n')
fprintf('Press Enter to proceed.\n')
pause
%</comment>

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

%% Simulation with perturbation on the target
opt_sim.dt = 0.02;
opt_sim.i_max = 3000;
opt_sim.tol = 0.1;
% Applying perturbations
opt_sim.perturbation.type='tcp'; %moving target continuously
opt_sim.perturbation.t0=1;
opt_sim.perturbation.tf=2;
opt_sim.perturbation.dx=[-20;20]; %we might as well randomize the perturbations

%<comment>
fprintf('\nNow we evaluate the model in the face of perturbation.\n')
fprintf('We start moving target from time t=1 sec till t=2 sec\n')
fprintf('with the constand velocity [20;-20] mm/s.\n')
fprintf('Press Enter to start simulator.\n')
pause
%</comment>

%running the simulator with the initial points of all demonstrations
fn_handle = @(x) GMR(Priors,Mu,Sigma,x,1:d,d+1:2*d);
[x_pert_tcp, xd_pert_tcp, tmp, xT_pert_tcp]=Simulation(x0_all,[],fn_handle,opt_sim);
plot_results(x_pert_tcp,xd_pert_tcp,Data,xT_pert_tcp,Mu,Sigma);

%% Simulation with perturbation on the target
opt_sim.dt = 0.02;
opt_sim.i_max = 3000;
opt_sim.tol = 0.1;
% Applying perturbations
opt_sim.perturbation.type='rdp'; %pushing robot instantly
opt_sim.perturbation.t0=1;
opt_sim.perturbation.dx=[-30;30]; %we might as well randomize the perturbations

%<comment>
fprintf('\nNow we evaluate the presence of discrete perturbation.\n')
fprintf('We pushed the robot''s end-effector at t=1 sec with the\n')
fprintf('displacement vector [-30;30] mm.\n')
fprintf('Press Enter to start simulator.\n')
pause
%</comment>

%running the simulator with the initial points of all demonstrations
fn_handle = @(x) GMR(Priors,Mu,Sigma,x,1:d,d+1:2*d);
[x_pert_rdp, xd_pert_rdp, tmp, xT_pert_rdp]=Simulation(x0_all,[],fn_handle,opt_sim);
plot_results(x_pert_rdp,xd_pert_rdp,Data,xT_pert_rdp,Mu,Sigma);






%=========================================================================%
function D = plot_results(x,xd,Data,xT,Mu,Sigma)
% plotting the result
figure('name','Results from Simulation','position',[265   200   520   720])
plot_legend; %plotting legend
sp(1)=subplot(3,1,1);
hold on; box on
plotGMM(Mu(1:2,:), Sigma(1:2,1:2,:), [0.6 1.0 0.6], 1,[0.6 1.0 0.6]);
plot(Data(1,:),Data(2,:),'r.')
xlabel('$\xi_1 (mm)$','interpreter','latex','fontsize',15);
ylabel('$\xi_2 (mm)$','interpreter','latex','fontsize',15);
title('Simulation Results')

sp(2)=subplot(3,1,2);
hold on; box on
plotGMM(Mu([1 3],:), Sigma([1 3],[1 3],:), [0.6 1.0 0.6], 1,[0.6 1.0 0.6]);
plot(Data(1,:),Data(3,:),'r.')
xlabel('$\xi_1 (mm)$','interpreter','latex','fontsize',15);
ylabel('$\dot{\xi}_1 (mm/s)$','interpreter','latex','fontsize',15);

sp(3)=subplot(3,1,3);
hold on; box on
plotGMM(Mu([2 4],:), Sigma([2 4],[2 4],:), [0.6 1.0 0.6], 1,[0.6 1.0 0.6]);
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
end
plot(sp(1),xT(1,end),xT(2,end),'k*','markersize',15,'linewidth',3)
plot(sp(2),xT(1,end),0,'k*','markersize',15,'linewidth',3)
plot(sp(3),xT(2,end),0,'k*','markersize',15,'linewidth',3)
D = axis(sp(1));


function [sl l]=plot_legend
sl=subplot(4,1,1);
hold on
plot(0,0,'k*','markersize',10,'linewidth',2);
plot([0 1],[0 1],'r:','linewidth',2);
plot(0,0,'kx','markersize',6,'linewidth',2);
% plot([0 1],[0 1],'k','linewidth',1);
patch([0 1],[0 1],[0.6 1 0.6],'linewidth',1);
plot([0 1],[0 1],'b','linewidth',1.5);
set(sl,'position',[1.2 0.76 0.1 0.1]);
l=legend('target','demonstrations','\mu','\Sigma','reproductions',...
        'orientation','horizontal','location','NorthOutside');
set(l,'position',[-0.0022    0.9642    1.0051    0.0347])

lc = get(l,'children');
r = diff(get(lc(4),'XData'))/2;
r = r(2);
center = mean([get(lc(4),'XData') get(lc(4),'YData')]);
t = linspace(-pi, pi, 40);
X = [cos(t);sin(t)] * r;
set(lc(4),'XData',X(1,:)+center(1),'YData',X(2,:)*7+center(2))
i_s=10;
x_tmp=get(lc(i_s),'XData');y_tmp=get(lc(i_s),'YData');
set(lc(i_s),'linestyle','none','marker','.','markersize',6,...
    'XData',linspace(x_tmp(1),x_tmp(2),8),'YData',linspace(y_tmp(1),y_tmp(2),8))