%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo Script for Drawing Data GMM-based LPV_DS Learning for paper:       %
%  'A Physically-Consistent Bayesian Non-Parametric Mixture Model for     %
%   Dynamical System Learning.'                                           %
% With this script you can load 2D toy trajectories or even real-world 
% trajectories acquired via kinesthetic taching and test the different    %
% GMM fitting approaches.                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Draw 2D Dataset with GUI  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc

% Create Figure
fig1 = figure('Color',[1 1 1]);
limits = [-6 0.5 -0.5 2];
axis(limits)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.55, 0.2646 0.4358]);
grid on

% Draw Reference Trajectories
[data, hp] = draw_mouse_data_on_DS(fig1, limits);

% Process Drawn Data for DS learning
[Data, Data_sh, att, x0_all, dt] = processDrawnData(data);

% Position/Velocity Trajectories
vel_samples = 10; vel_size = 0.5; 
[h_data, h_att, h_vel] = plot_reference_trajectories_DS(Data, att, vel_samples, vel_size);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% To use this dataset go to the demo_loadData_*.m scripts %
%% and start with the[Step 2] block of code                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%