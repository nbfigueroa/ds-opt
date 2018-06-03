function Tutorial_Obstacle_Avoidance
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Copyright (c) 2011 S. Mohammad Khansari-Zadeh, LASA Lab, EPFL,   %%%
%%%          CH-1015 Lausanne, Switzerland, http://lasa.epfl.ch         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The program is submitted as an accompanying material to this submission:
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
disp('In this tutorial, we show examples of using the proposed obstacle')
disp('avoidance algorithm with different dynamical systems (DS). Throughout')
disp('this demo, both the DS and the obstacle equations are given. The file')
disp('''obs_modulation_ellipsoid.m'' in the ''lib_obstacle_avoidance'' folder,')
disp('is the only file that is necessary to perform the obstacle avoidance.')
disp('All other files are just needed for the illustrative purpose.')
disp('Please FIRST run ''Tutotial_Obstacle_Avoidance.m'' tutorial to get familiar.')
disp('with the format with which we define obstacles.')
disp('press any key to continue ...')
pause

%% first demo
disp(' ')
disp(' ')
disp('In the first demo, we consider a globally asymptotically stable DS,')
disp('which is defined by:')
disp('  xd = - x(1,:);')
disp('  xd(2,:) = - x(1,:) .* cos(x(1,:)) - x(2,:);')
disp('This DS has a unique attractor at the origin.')
disp('press any key to draw the streamlines of this DS ...')
pause
fn_handle = @(x) globally_stable_DS(x); %defining the function handle
x0 = [-15*ones(1,15);linspace(-5,5,15)]; %set of initial points
% A set of parameters that should be defined for the simulation
opt_sim.dt = 0.02; %integration time steps
opt_sim.i_max = 1000; %maximum number of iterations
opt_sim.tol = 0.05; %convergence tolerance
opt_sim.plot = true; %enabling the animation
opt_sim.model = 1; %first order ordinary differential equation
opt_sim.obstacle = []; %no obstacle is defined
fig(1) = figure('name','First demo: Streamlines of the original DS','position',[100 550 560 420]);
opt_sim.figure = fig(1);
Simulation(x0,[],fn_handle,opt_sim);
disp(' ')
disp(' ')
disp('Now, let us evaluate the system behavior in the presence of a single')
disp('obstacle, that is defind by:')
disp('obs{1}.a  = [1.2 1.2')
disp('             0.4 1.0]')
disp('obs{1}.p  = [2 1')
disp('             1 1]')
disp('obs{1}.x0 = [-8;0]')
disp('obs{1}.sf = [1.2;1.2]')
disp('obs{1}.th_r = 0')
disp('obs{1}.partition = [-pi 0')
disp('                      0 pi]')
disp('press any key to draw the streamlines in the presence of the obstacle ...')
pause
clear obs;
obs{1}.a = [1.2 1.2;0.4 1];
obs{1}.p = [2 1;1 1];
obs{1}.partition = [-pi 0;0 pi];
obs{1}.x0 = [-8;0];
obs{1}.sf = [1.2;1.2];
obs{1}.th_r = 0*pi/180;
opt_sim.obstacle = obs;
fig(2) = figure('name','First demo: Streamlines of the modulated DS','position',[660 550 560 420]);
opt_sim.figure = fig(2);
Simulation(x0,[],fn_handle,opt_sim);
disp('press any key to continue ...')
pause

% adding more obstacles
disp(' ')
disp(' ')
disp('Now let us add two more obstacles to the previous example. We use the')
disp('same geometry for the new obstacles, but we place them in different')
disp('positions, and also rotate them with different angles.')
disp('The second obstacle:')
disp('  obs{2} = obs{1}')
disp('  obs{2}.x0 = [-12;3]')
disp('  obs{2}.th_r = 90*pi/180')
disp('The third obstacle:')
disp('  obs{3} = obs{1}')
disp('  obs{3}.x0 = [-12;-3]')
disp('  obs{3}.th_r = -90*pi/180')
obs{2} = obs{1};
obs{2}.x0 = [-12;3];
obs{2}.th_r = 90*pi/180;
obs{3} = obs{1};
obs{3}.x0 = [-12;-3];
obs{3}.th_r = -90*pi/180;
disp('press any key to draw the streamlines in the presence of three obstacles ...')
pause
opt_sim.obstacle = obs;
fig(3) = figure('name','First demo: Multiple obstacle avoidance','position',[380 50 560 420]);
opt_sim.figure = fig(3);
Simulation(x0,[],fn_handle,opt_sim);
disp('press any key to continue ...')
pause
close(fig)
%% second demo
disp(' ')
disp(' ')
disp('In the second demo, we consider a stable limit cycle that drive the robot')
disp('motion:')
disp('  xd = x(2,:);')
disp('  xd(2,:) = -x(1,:) + 0.9*(1 - x(1,:).^2).*x(2,:);')
disp('press any key to draw the streamlines of this DS ...')
pause
fn_handle = @(x) limitcycle_DS(x); %defining the function handle
x0 = [zeros(1,10);linspace(0,3,10)]; %set of initial points
% A set of parameters that should be defined for the simulation
opt_sim.dt = 0.01; %integration time steps
opt_sim.i_max = 1000; %maximum number of iterations
opt_sim.tol = 0.05; %convergence tolerance
opt_sim.plot = true; %enabling the animation
opt_sim.obstacle = []; %no obstacle is defined
fig=[];
fig(1) = figure('name','Second demo: Streamlines of the original DS','position',[100 550 560 420]);
opt_sim.figure = fig(1);
Simulation(x0,[],fn_handle,opt_sim);

disp(' ')
disp(' ')
disp('For simplicity, we use the same geometry as the previous demo, but')
disp('with different positions. In this example, we consider two obstacles.')
disp('press any key to draw the streamlines in the presence of the obstacles ...')
pause
obs(3) = []; %deleting the third obstacle
obs{1}.x0 = [-3;-1];
obs{2}.x0 = [3;1];

opt_sim.obstacle = obs;
fig(2) = figure('name','Second demo: Multiple obstacle avoidance','position',[660 550 560 420]);
opt_sim.figure = fig(2);
Simulation(x0,[],fn_handle,opt_sim);
disp('press any key to continue ...')
pause
close(fig)
%% third demo
disp(' ')
disp(' ')
disp('In the third demo, we go through how to modify the characteristics')
disp('of the avoidance through the parameters:')
disp('   - reactivity    (.rho)')
disp('   - Tail Effect   (.tailEffect)')
disp('To have a better visualisation, we consider a simple dynamical systems:')
disp('  xd = [-x(1,:);0];')
disp('press any key to draw the streamlines of this DS ...')
pause
fn_handle = @(x) [-x(1,:);zeros(1,size(x,2))]; %defining the function handle
x0 = [-10*ones(1,20);linspace(-3,3,20)]; %set of initial points
% A set of parameters that should be defined for the simulation
opt_sim.dt = 0.025; %integration time steps
opt_sim.i_max = 1000; %maximum number of iterations
opt_sim.tol = 0.05; %convergence tolerance
opt_sim.plot = true; %enabling the animation
opt_sim.obstacle = []; %no obstacle is defined
fig = [];
fig(1) = figure('name','Third demo: Streamlines of the original DS','position',[100 550 560 420]);
opt_sim.figure = fig(1);
Simulation(x0,[],fn_handle,opt_sim);

disp(' ')
disp(' ')
disp('For simplicity, we use the same geometry as the previous demo, but')
disp('with different positions. In this example, we consider one obstacle.')
disp('press any key to draw the streamlines in the presence of the obstacle ...')
pause
obs(2) = []; %deleting the third obstacle
obs{1}.x0 = [-5;-1];
opt_sim.obstacle = obs;
fig(2) = figure('name','Third demo: Streamlines of the modulated DS','position',[660 550 560 420]);
opt_sim.figure = fig(2);
Simulation(x0,[],fn_handle,opt_sim);
disp('press any key to continue ...')
pause

% Changing the reactivity
disp(' ')
disp(' ')
disp('Now let us change the obstacle property so that the motion reacts')
disp('earlier to the presence of the robot. For this, we increase the')
disp('reactivity factor from 1 to 3. We add the following line:')
disp('   obs{1}.rho = 3')
obs{1}.rho = 3;
obs{1}.tailEffect = true;
disp('press any key to draw the streamlines ...')
pause
opt_sim.obstacle = obs;
fig(3) = figure('name','Third demo: Effect of the reactivity parameter','position',[100 50 560 420]);
opt_sim.figure = fig(3);
Simulation(x0,[],fn_handle,opt_sim);
disp('press any key to continue ...')
pause

% Removing the taileffect
disp(' ')
disp(' ')
disp('Now let us change another obstacle property. As you saw in previous')
disp('examples, the obstacle modifies the motion even when the it is moving')
disp('moving away from the obstacle. If you don''t like this effect, you just')
disp('need to tell it to the module through the following line:')
disp('   obs{1}.tailEffect = false')
disp('By setting the tailEffect property to the false, the motion is only')
disp('modified when it approaches the obstacle.')
obs{1}.tailEffect = false;
disp('press any key to draw the streamlines ...')
pause
opt_sim.obstacle = obs;
fig(4) = figure('name','Third demo: Removing the tail effect','position',[660 50 560 420]);
opt_sim.figure = fig(4);
Simulation(x0,[],fn_handle,opt_sim);
disp('press any key to exit.')
pause
close(fig);
disp('The tutorial ended.')