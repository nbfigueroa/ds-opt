% This script requires the SEDS package to be on MATLAB path.
%
% This is an example script demonstrating the process to learn the coupling 
% model in CDS. Once the coupling is learned, it can be combined with the
% individual dynamic (SEDS) models to generate the trajectories of the overall
% system in a coupled manner.
% There are two examples attached. First is a 2-D example where the first dimension
% is master and the other slave. Second is a 3-D example where the first
% dimension is master and the latter 2 dimensions are slave. 
%
% One could also create new 2-D data using the mouse.
%
% Format of the data is the same as SEDS data format. Type 'doc preprocess_demos'
% in the MATLAB command window for more information.

%% Example data sets
close all

%Ex. 1
load('data/ex1_Data.mat', 'full_data');                    
master_dim = 1;                 % indices of the master variables
slave_dim = 2;                  % indices of the slave variables
coupl_func = @(x) abs(x);       % coupling function

% %Ex. 2
% load ('data/ex2_Data.mat','full_data');                    
% master_dim = [3];
% slave_dim = 1:2;
% coupl_func = @(x) norm(x);

% Preprocess data ==> centering, smoothing, truncating
[x0 , xT, Data, index] = preprocess_demos(full_data,0.05,0.000001);

% Plotting raw data for visualization
figure(1)
clf
hold on
if(size(Data,1)/2 == 2)
    plot(Data(1,:)', Data(2,:)','k.');
    xlabel('$x_1$','interpreter','latex','Fontsize',25);
    ylabel('$x_2$','interpreter','latex','Fontsize',25);
    grid on
else
    plot3(Data(1,:)', Data(2,:)', Data(3,:)','k.');
    xlabel('$x_1$','interpreter','latex','Fontsize',25);
    ylabel('$x_2$','interpreter','latex','Fontsize',25);
    zlabel('$x_3$','interpreter','latex','Fontsize',25);
    grid on
    view(3);
end
hl=legend('Data');
set(hl,'interpreter','latex','Fontsize',20);
axis equal
axis auto

% Learn coupling model
K = 5;                  % Number of gaussians for learning
structGMM = learn_coupling(Data, master_dim, slave_dim, coupl_func, K);

save('data/coupling_model.mat','structGMM');

%% Plotting Model

% Applying coupling function on the master states
master_data=[];
for i=1:size(Data,2)
master_data = [master_data, coupl_func(Data(master_dim,i))];
end

for i=1:size(master_data,1)   % If output of the coupling function is multi-dimensional
    figure(i)
    clf
    for j=1:length(slave_dim) % All slave dimensions vs master_i in figure i
        subplot(1,length(slave_dim),j);
        cla; hold on; box on
        plotGMM(structGMM.Mu([i,j+size(master_data,1)],:), structGMM.Sigma([i,j+size(master_data,1)],...
            [i,j+size(master_data,1)],:), [0.6 1.0 0.6], 1,[0.6 1.0 0.6]);
        plot(master_data(i,:)',Data(slave_dim(j),:)','r.');
        x_in = linspace(0,max(master_data(i,:)),100);
        [tmp, sig_tmp] = GMR(structGMM.Priors, structGMM.Mu, structGMM.Sigma, x_in,i,j+size(master_data,1));
        plot(x_in', tmp','k-','Linewidth',2);
        xlabel(['$M_' int2str(i) '$'], 'Fontsize',25, 'interpreter','latex');
        ylabel(['$S_' int2str(j) '$'], 'Fontsize',25, 'interpreter','latex');
    end
end

%% Coupled dynamics/trajectories in the original space

% Load dynamic models of master and slave systems

%Ex.1 
load data/ex1_masterDyn
m_gmm = structGMM;
load data/ex1_slaveDyn
s_gmm = structGMM;
load data/coupling_model
cpl_gmm = structGMM;
cpl_func = structGMM.cplfunc;

% %Ex. 2
% load data/ex2_masterDyn.mat
% m_gmm = structGMM;
% load data/ex2_slaveDyn
% s_gmm = structGMM;
% load data/coupling_model
% cpl_gmm = structGMM;
% cpl_func = structGMM.cplfunc;


% Checking for compatible dimensionalities of different models
if(size(m_gmm.Mu,1)/2 ~= length(cpl_func(zeros(length(cpl_gmm.master_dim),1))))
    disp('Bad dimensions of the master dynamic model!');
end
if(size(s_gmm.Mu,1)/2 ~= length(cpl_gmm.slave_dim))
    disp('Bad dimensions of the slave dynamic model!');
end

% Dimensionality of the state space
dim = (size(m_gmm.Mu,1)+size(s_gmm.Mu,1))/2;
if(dim==2)
    %Open parameters ==> Reaction speed and amplitude
    alpha=[1.0;5.0]; beta=[1.0;5.0];
    [tx,ty] = meshgrid(-1:0.2:2, -2:0.2:2);
    vx1=[];vy1=[];
    vx2=[];vy2=[];
    vx3=[];vy3=[];
    
    for i=1:size(tx,1)
        for j=1:size(tx,2)
            point = [tx(i,j);ty(i,j)];
            tmp = getCDSVelocity(m_gmm, s_gmm, cpl_gmm, point, alpha(1), beta(1));
            vx1(i,j)=tmp(1);
            vy1(i,j)=tmp(2);
            
            tmp = getCDSVelocity(m_gmm, s_gmm, cpl_gmm, point, alpha(2), beta(1));
            vx2(i,j)=tmp(1);
            vy2(i,j)=tmp(2);
            
            tmp = getCDSVelocity(m_gmm, s_gmm, cpl_gmm, point, alpha(2), beta(2));
            vx3(i,j)=tmp(1);
            vy3(i,j)=tmp(2);
        end
    end
    % Plotting streamlines for the above parameter values
    figure(2)
    clf
    subplot(1,3,1);
    cla
    hold on
    box on
    streamslice(tx,ty,vx1,vy1,1.5);
    plot(Data(1,:)',Data(2,:)','k.');
    plot(0,0,'ko','Markerfacecolor','k', 'Markersize',12);
    title('$\alpha=1.0;\beta=1.0$','Fontsize',25,'interpreter','latex');
    xlabel('$x_1$','interpreter','latex','Fontsize',25);
    ylabel('$x_2$','interpreter','latex','Fontsize',25);
    axis equal
    axis tight
    
    subplot(1,3,2);
    cla
    hold on
    box on
    streamslice(tx,ty,vx2,vy2,1.5);
    plot(Data(1,:)',Data(2,:)','k.');
    plot(0,0,'ko','Markerfacecolor','k', 'Markersize',12);
    title('$\alpha=5.0;\beta=1.0$','Fontsize',25,'interpreter','latex');
    xlabel('$x_1$','interpreter','latex','Fontsize',25);
    ylabel('$x_2$','interpreter','latex','Fontsize',25);
    axis equal
    axis tight
    
    subplot(1,3,3);
    cla
    hold on
    box on
    ht=streamslice(tx,ty,vx3,vy3,1.5);
    hd=plot(Data(1,:)',Data(2,:)','k.');
    title('$\alpha=5.0;\beta=5.0$','Fontsize',25,'interpreter','latex');
    xlabel('$x_1$','interpreter','latex','Fontsize',25);
    ylabel('$x_2$','interpreter','latex','Fontsize',25);
    
    ht2=plot(0,0,'ko','Markerfacecolor','k', 'Markersize',12);
    hl=legend([hd, ht2, ht(1)],'Data','Target','Trajectories');
    set(hl, 'interpreter','latex','Fontsize',20);
    axis equal
    axis tight
    
else if(dim==3)
        %Open parameters ==> Reaction speed and amplitude
        alpha=[0.5;1.0;2.0];beta=[0.5;1.0;1.0];
        figure(2)
        cla
        hold on
        grid on
        box on
        title('Resulting trajectories with coupled dynamics','interpreter','latex','Fontsize',20);
        ht1=plot3(Data(1,:)', Data(2,:)', Data(3,:)','k.');
        axis equal
        
        % Initial points for simulating trajecotries
        x=[rand(2,1)*0.8 ;Data(3,1)/6];
        x2=[0.5+rand(2,1);Data(3,1)/2];
        x3=[0.5-rand(2,1);Data(3,1)/1];
        
        v=getCDSVelocity(m_gmm, s_gmm, cpl_gmm, x(:,end), alpha(1), beta(1));
        v2=getCDSVelocity(m_gmm, s_gmm, cpl_gmm, x2(:,end), alpha(2), beta(2));
        v3=getCDSVelocity(m_gmm, s_gmm, cpl_gmm, x3(:,end), alpha(3), beta(3));
        
        plot3(x(1,1),x(2,1),x(3,1),'ko','Linewidth',3, 'Markersize',12);
        plot3(x2(1,1),x2(2,1),x2(3,1),'ko','Linewidth',3, 'Markersize',12);
        ht2=plot3(x3(1,1),x3(2,1),x3(3,1),'ko','Linewidth',3, 'Markersize',12);
        
        h=plot3(x(1,1),x(2,1),x(3,1));
        h2=plot3(x2(1,1),x2(2,1),x2(3,1));
        h3=plot3(x3(1,1),x3(2,1),x3(3,1));
        
        xlabel('$x_1$','interpreter','latex','Fontsize',25);
        ylabel('$x_2$','interpreter','latex','Fontsize',25);
        zlabel('$x_3$','interpreter','latex','Fontsize',25);
        % Target point
        ht3=plot3(0,0,0,'ko','Markerfacecolor','k', 'Markersize',12);
        hl=legend([ht1,ht2,ht3,h3],'Data','Starting points','Target','Trajectories');
        set(hl, 'interpreter','latex','Fontsize',20);
        axis auto
        
        % Simulating trajectories starting from various points
        while(norm(v) > 1e-2 || norm(v2) > 1e-2 || norm(v3) > 1e-2)
            v = getCDSVelocity(m_gmm, s_gmm, cpl_gmm, x(:,end), alpha(1), beta(1));
            x(:,end+1)=x(:,end)+v*0.05;
            set(h, 'XData',x(1,:)', 'YData',x(2,:)', 'ZData',x(3,:)');
            
            v2 = getCDSVelocity(m_gmm, s_gmm, cpl_gmm, x2(:,end), alpha(2), beta(2));
            x2(:,end+1)=x2(:,end)+v2*0.05;
            set(h2, 'XData',x2(1,:)', 'YData',x2(2,:)', 'ZData',x2(3,:)');
            
            v3 = getCDSVelocity(m_gmm, s_gmm, cpl_gmm, x3(:,end), alpha(3), beta(3));
            x3(:,end+1)=x3(:,end)+v3*0.05;
            set(h3, 'XData',x3(1,:)', 'YData',x3(2,:)', 'ZData',x3(3,:)');
            pause(0.01);
        end
    else        % For higher dimensional models. Can be also used for 2-D/3-D models.
        
        %Set parameters
        alpha = 2.0; beta=1.0;
        
        figure(2)
        x=[];
        h=[];
        v=[];
        for i=1:dim
            subplot(dim,1,i);
            cla
            hold on
            grid on
            for j=1:length(index)-1
                plot(Data(i,index(j):index(j+1)-1)','k.');
            end
            x(i,1) = Data(i,1);
            plot(x(i,1),'ko','Linewidth',3, 'Markersize',12);
            h(i) = plot(x(i));
            ylabel(['$x_' int2str(i) '$'], 'interpreter','latex','Fontsize',25);
            xlabel('Timestep', 'interpreter','latex','Fontsize',25);
        end
        
        v = getCDSVelocity(m_gmm, s_gmm, cpl_gmm, x(:,end), alpha, beta);
        while(norm(v) > 1e-2)
            v = getCDSVelocity(m_gmm, s_gmm, cpl_gmm, x(:,end), alpha, beta);
            x(:,end+1)=x(:,end)+v*0.05;
            for i=1:dim
                set(h(i), 'YData',x(i,:)');
            end
            pause(0.01);
        end
    end
end

