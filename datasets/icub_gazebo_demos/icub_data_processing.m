%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Data Processing Script      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all
pkg_dir         = '/home/nbfigueroa/Dropbox/PhD_papers/CoRL-2018/code/ds-opt/';
load(strcat(pkg_dir,'datasets/icub_gazebo_demos/raw_data'))

data = []; 
odata = [];
window_size = 151; crop_size = (window_size+1)/2; 
dt = mean(abs(diff(raw_data{1}(1,:))));

for i=1:length(raw_data)
        dx_nth = sgolay_time_derivatives(raw_data{i}(2:3,:)', dt, 2, 2, window_size);
        X     = dx_nth(:,:,1)';
        X_dot = dx_nth(:,:,2)';               
        data{i} = [X; X_dot];
        
        % Extract orientation data too
        odata{i} = raw_data{i}(4,crop_size:end-crop_size);
end

% Trajectories to use
left_traj = 1;
if ~left_traj
    data(4:end) = [];
    odata(4:end) = [];
end

% To use sub-sample measurments
sub_sample = 5;
[Data, Data_sh, att, x0_all, dt, data] = processDataStructure2(data, sub_sample);

% Position/Velocity Trajectories
vel_samples = 80; vel_size = 0.75; 
[h_data, h_att, h_vel] = plot_reference_trajectories_DS(Data, att, vel_samples, vel_size);
axis equal

% Draw Obstacles
rectangle('Position',[-1 1 6 1], 'FaceColor',[.85 .85 .85]); hold on;
h_att = scatter(att(1),att(2), 150, [0 0 0],'d','Linewidth',2); hold on;



