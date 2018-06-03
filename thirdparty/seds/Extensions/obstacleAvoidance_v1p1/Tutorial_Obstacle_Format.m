function Tutorial_Obstacle_Format
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Copyright (c) 2011 S. Mohammad Khansari-Zadeh, LASA Lab, EPFL,   %%%
%%%          CH-1015 Lausanne, Switzerland, http://lasa.epfl.ch         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is submitted as an accompanying material to this submission:
% 
%     S.M. Khansari Zadeh and A. Billard, "A Dynamical System Approach to
%     Realtime Obstacle Avoidance", Autonomous Robots, 2012
%
%% preparing the obstacle avoidance module
%adding the obstacle avoidance folder to the MATLAB path directories
if isempty(regexp(path,['lib_obstacle_avoidance' pathsep], 'once'))
    addpath([pwd, '/lib_obstacle_avoidance']);
end
%adding the example folder to the MATLAB path directories
if isempty(regexp(path,['Example_DSs' pathsep], 'once'))
    addpath([pwd, '/Example_DSs']);
end
%% initial comments
clc
disp('The aim of this tutorial is to help readers to properly define the')
disp('analytical model of obstacles in the format that is consistent with')
disp('the provided obstacle avoidance library.')
disp('press any key to continue ...')
pause

disp(' ')
disp(' ')
disp('To avoid struggling with the MATLAB symbolic toolbox, and for the sake')
disp('of simplicity, the provided source code can only be used to avoid obstacles')
disp('in the form of:')
disp('          \Gamma(xt):   \sum_{i=1}^d (xt_i/a_i)^(2p_i) = 1')
disp('For other forms of obstacle shape, it is necessary to modify the files')
disp('in the ''lib_obstacle_avoidance'' folder to properly compute the matrix')
disp('of eigenvalues and eigenvectors of the dynamic modulation matrix.')
disp('press any key to continue ...')
pause
%% obstacle format
disp(' ')
disp(' ')
disp('The provided obstacle avoidance library accepts the obstacle data in a')
disp('special form. Each obstacle should be defined with a structure, and the')
disp('properties of all the obstacles need to be merged into a cell variable.')
disp('For example, for a cell array ''obs'', obs{i} provides the properties of')
disp('the i-th obstacle. Each obs{i} is a structure variable, and should contains')
disp('at least the following properties:')
disp('obs{i}.a  : The ellipsoid length scale (e.g. [1;1] for a 2D circle)')
disp('obs{i}.p  : The ellipsoid power terms (e.g. [1;1] for a 2D circle)')
disp('obs{i}.x0 : The ellipsoid reference point (e.g. [0;0])')
disp('obs{i}.sf (optional, default = 1): The safety factor (e.g. [1.5;1.5]).')
disp('obs{i}.tailEffect (optional, default = true): If it is set to true, the obstacle');
disp('            modifies the motion even when the it is moving away from the obstacle.');
disp('obs{i}.rho (optional, default = 1.0): Sets the reactivity parameter. The larger');
disp('            the rho, the earlier the robot responds to the presence of the obstacle.');
disp('obs{i}.th_r (optional, default = 0.0): The obstacle rotation with respect to ');
disp('            the global frame of reference (e.g. 0 rad)')
disp('obs{i}.partition (optional, default = [-pi pi]): Defines the partition of the ellipsoid')
disp('We will elaborate more the ''.partition'' property later on.')
disp('press any key to continue ...')
pause
%% example of an ellipsoid obstacle
disp(' ')
disp(' ')
disp('For example these are the properties to define an ellipse with major ')
disp('axes 1 and 2, centered at the origin, with the safety factor of 1.2:')
disp('obs{1}.a  = [1;2]')
disp('obs{1}.p  = [1;1]')
disp('obs{1}.x0 = [0;0]')
disp('obs{1}.sf = [1.2;1.2]')
disp('obs{1}.th_r = 0*pi/180')
disp('obs{1}.partition = [-pi pi]')
disp('Recall, the obstacle formulation is in the form of: \sum_{i=1}^d (xt_i/a_i)^(2p_i) = 1')
disp('press any key to draw the obstacle ...')
pause

clear obs;
obs{1}.a  = [1;2];
obs{1}.p  = [1;1];
obs{1}.x0 = [0;0];
obs{1}.sf = [1.2;1.2];
obs{1}.th_r = 0*pi/180;
obs{1}.partition = [-pi pi];

np = 200; %np is the number of points used to draw the obstacle
[x_obs x_obs_sf] = obs_draw_ellipsoid(obs,np); 
figure('position',[300 300 560 420])
box on;grid on;hold on;axis equal
for n=1:length(obs)
    patch(x_obs(1,:,n),x_obs(2,:,n),[0.6 1 0.6])
    plot(x_obs_sf(1,:,n),x_obs_sf(2,:,n),'k--','linewidth',0.5)
end
axis([-6 6 -4 4])
disp('press any key to continue ...')
pause
%% example on safety factor
disp(' ')
disp(' ')
disp('By changing .sf property, we can have different safety factor')
disp('along each axis:')
disp('obs{1}.sf = [2;1.2]')
disp('press any key to draw the obstacle ...')
pause

obs{1}.sf = [2;1.2];

close(gcf)
np = 200; %np is the number of points used to draw the obstacle
x_obs = obs_draw_ellipsoid(obs,np); 
figure('position',[300 300 560 420])
box on;grid on;hold on;axis equal
[x_obs x_obs_sf] = obs_draw_ellipsoid(obs,np); 
for n=1:length(obs)
    patch(x_obs(1,:,n),x_obs(2,:,n),[0.6 1 0.6])
    plot(x_obs_sf(1,:,n),x_obs_sf(2,:,n),'k--','linewidth',0.5)
end
axis([-6 6 -4 4])
disp('press any key to continue ...')
pause
%% example on defining a squarish obstacle
disp(' ')
disp(' ')
disp('By changing .p property, we can also define a more squarish obstacle:')
disp('obs{1}.p = [3;3]')
disp('press any key to draw the obstacle ...')
pause

obs{1}.p = [3;3];

close(gcf)
np = 200; %np is the number of points used to draw the obstacle
[x_obs x_obs_sf] = obs_draw_ellipsoid(obs,np); 
figure('position',[300 300 560 420])
box on;grid on;hold on;axis equal
for n=1:length(obs)
    patch(x_obs(1,:,n),x_obs(2,:,n),[0.6 1 0.6])
    plot(x_obs_sf(1,:,n),x_obs_sf(2,:,n),'k--','linewidth',0.5)
end
axis([-6 6 -4 4])
disp('press any key to continue ...')
pause
%% example on rotating the obstacle
disp(' ')
disp(' ')
disp('We can also rotate the obstacle by changing .th_r property:')
disp('obs{1}.th_r = 60*pi/180')
disp('press any key to draw the obstacle ...')
pause

obs{1}.th_r = 60*pi/180;

close(gcf)
np = 200; %np is the number of points used to draw the obstacle
[x_obs x_obs_sf] = obs_draw_ellipsoid(obs,np); 
figure('position',[300 300 560 420])
box on;grid on;hold on;axis equal
for n=1:length(obs)
    patch(x_obs(1,:,n),x_obs(2,:,n),[0.6 1 0.6])
    plot(x_obs_sf(1,:,n),x_obs_sf(2,:,n),'k--','linewidth',0.5)
end
axis([-6 6 -4 4])
disp('press any key to continue ...')
pause
%% example one defining a complex obstacle via partitioning
disp(' ')
disp(' ')
disp('It is possible to generate more complex shape obstacles, by partitioning the space.')
disp('For example, suppose that we want to define an obstacle with:')
disp('  x^2 + y^2 = 1         when y<0')
disp('  x^2 + (3*y)^4 = 1     when y>=0')
disp('The can model the above formulation by using:')
disp('obs{1}.a  = [1 1')
disp('             1 3]')
disp('obs{1}.p  = [1 1')
disp('             1 2]')
disp('obs{1}.partition = [-pi 0')
disp('                     0 pi]')
disp('The number of rows in .partition indicates the number of partitions that is used to')
disp('define the obstacle. The i-th column in .a and .p provides the obstacle property for')
disp('the i-th partition defined in the i-th row of .partition. Note that the number of')
disp('columns in .a and .p should be equal to the number of rows in .partition.')
disp('The user should ensure the continuity and the differentiability of the ')
disp('formulation at the transition points across the partitions.')
disp('press any key to draw the new obstacle ...')
pause
obs{1}.a  = [1 1;1 3];
obs{1}.p  = [1 1;1 2];
obs{1}.partition = [-pi 0;0 pi];

close(gcf)
np = 200; %np is the number of points used to draw the obstacle
[x_obs x_obs_sf] = obs_draw_ellipsoid(obs,np); 
figure('position',[300 300 560 420])
box on;grid on;hold on;axis equal
for n=1:length(obs)
    patch(x_obs(1,:,n),x_obs(2,:,n),[0.6 1 0.6])
    plot(x_obs_sf(1,:,n),x_obs_sf(2,:,n),'k--','linewidth',0.5)
end
axis([-6 6 -4 4])
disp('press any key to continue ...')
pause
%% example on defining multiple obstacles
disp(' ')
disp(' ')
disp('Extension to multiple obstacles can be simply obtained by defining obs{2}, obs{3}, ...')
disp('For example, now we just add a circular obstacle next to the obstacle we have already defined.')
disp('The new obstacle property is:')
disp('obs{2}.a  = [1.5;1.5]')
disp('obs{2}.p  = [1;1]')
disp('obs{2}.x0 = [4;-2]')
disp('obs{2}.sf = [1.2;1.2]')
disp('obs{2}.th_r = 0*pi/180')
disp('obs{2}.partition = [-pi pi]')
disp('press any key to draw the obstacles ...')
pause

obs{2}.a  = [1.5;1.5];
obs{2}.p  = [1;1];
obs{2}.x0 = [4;-2];
obs{2}.sf = [1.2;1.2];
obs{2}.th_r = 0*pi/180;
obs{2}.partition = [-pi pi];

np = 200; %np is the number of points used to draw the obstacle
[x_obs x_obs_sf] = obs_draw_ellipsoid(obs,np); 
close(gcf)
figure('position',[300 300 560 420])
box on;grid on;hold on;axis equal
for n=1:length(obs)
    patch(x_obs(1,:,n),x_obs(2,:,n),[0.6 1 0.6])
    plot(x_obs_sf(1,:,n),x_obs_sf(2,:,n),'k--','linewidth',0.5)
end
axis([-6 6 -4 4])
disp('press any key exit.')
pause
%% example on defining a 3D obstacle
disp(' ')
disp(' ')
disp('Defining multi-dimensional obstacles is straight forward. For example, ')
disp('the following properties can be used to define a 3D ellipsoid obstacle:')
disp('obs{2}.a  = [1;2;3]')
disp('obs{2}.p  = [3;2;2]')
disp('obs{2}.x0 = [0;0;0]')
disp('obs{2}.sf = [1;1;1]')
disp('obs{2}.th_r = [0 0 0]*pi/180')
disp('obs{2}.partition = [-pi pi]')
disp('press any key to draw the obstacles ...')
pause
clear obs
obs{1}.a  = [1;2;3];
obs{1}.p  = [3;2;2];
obs{1}.x0 = [0;0;0];
obs{1}.sf = [1;1;1];
obs{1}.th_r = [0 0 0]*pi/180;
obs{1}.partition = [-pi pi];

n_theta = 20;
n_phi = 15;
[x_obs x_obs_sf] = obs_draw_ellipsoid(obs,[n_theta n_phi]);
close(gcf)
figure('position',[300 300 560 420])
grid on;hold on;axis equal
for n=1:length(obs)
    h = surf(reshape(x_obs(1,:,n),n_phi,n_theta), reshape(x_obs(2,:,n),n_phi,n_theta), reshape(x_obs(3,:,n),n_phi,n_theta));
    set(h,'FaceColor',[0.6 1 0.6],'linewidth',0.1)
end
view(-50,38)
axis([-6 6 -4 4 -4 4])
%%
disp(' ')
disp(' ')
disp('However, as you have noticed, there are a few differences when defining')
disp('obstacles in higher dimensions. For example, for 3D models, .th_r property')
disp('has three components. Each component corresponds to one of the Euler angles')
disp('that is requires to rotate an object in 3D space. The first, second, and')
disp('third component of .th_r are respectively rotation around the x, y, and z')
disp('axes respectively.')
disp('For example, let us modify obs{1}.th_r to [-30 -45 0]*pi/180.')
obs{1}.th_r = [-30 -45 0]*pi/180;
disp('press any key to draw the rotated obstacle ...')
pause

[x_obs x_obs_sf] = obs_draw_ellipsoid(obs,[n_theta n_phi]);
close(gcf)
figure('position',[300 300 560 420])
grid on;hold on;axis equal
for n=1:length(obs)
    h = surf(reshape(x_obs(1,:,n),n_phi,n_theta), reshape(x_obs(2,:,n),n_phi,n_theta), reshape(x_obs(3,:,n),n_phi,n_theta));
    set(h,'FaceColor',[0.6 1 0.6],'linewidth',0.1)
end
view(-50,38)
axis([-6 6 -4 4 -4 4])
disp('press any to continue.')
pause
disp(' ')
disp(' ')
disp('Please note that the obstacle rotation is not supported in dimensions')
disp('higher than 3D.')
disp('press any key to continue ...')
pause
%%
disp(' ')
disp(' ')
disp('And finally here is a 3D example defined by using two partitions:')
disp('obs{1}.a  = [1.5 1.5')
disp('             2.1 1.5')
disp('             0.8 0.8]')
disp('obs{1}.p  = [1 1')
disp('             1 2')
disp('             1 1]')
disp('obs{1}.x0 = [0;0;0]')
disp('obs{1}.sf = [1;1;1]')
disp('obs{1}.th_r = [0 0 0]*pi/180')
disp('obs{1}.partition = [-pi 0')
disp('                      0 pi]')
obs{1}.a  = [1.5 1.5
             2.1 1.5
             0.8 0.8];
obs{1}.p  = [1 1
             1 2
             1 1];
obs{1}.x0 = [0;0;0];
obs{1}.sf = [1;1;1];
obs{1}.th_r = [0 0 0]*pi/180;
obs{1}.partition = [-pi 0
                      0 pi];

disp('press any key to draw the obstacle ...')
pause

[x_obs x_obs_sf] = obs_draw_ellipsoid(obs,[n_theta n_phi]);
close(gcf)
figure('position',[300 300 560 420])
grid on;hold on;axis equal
for n=1:length(obs)
    h = surf(reshape(x_obs(1,:,n),n_phi,n_theta), reshape(x_obs(2,:,n),n_phi,n_theta), reshape(x_obs(3,:,n),n_phi,n_theta));
    set(h,'FaceColor',[0.6 1 0.6],'linewidth',0.1)
end
view(-50,38)
axis([-6 6 -4 4 -4 4])
disp('press any key exit.')
pause
close(gcf)
disp('The tutorial ended.')