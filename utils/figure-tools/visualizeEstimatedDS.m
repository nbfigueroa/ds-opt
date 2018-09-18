function [hd, hs, hr, x_sim] = visualizeEstimatedDS(Xi_ref, ds_fun, ds_plot_options)
fig1 = figure('Color',[1 1 1]);
M = size(Xi_ref,1);

% Parse Options
plot_repr   = ds_plot_options.sim_traj;
x0_all      = ds_plot_options.x0_all;
init_type   = ds_plot_options.init_type;
nb_pnts     = ds_plot_options.nb_points;
plot_volume = ds_plot_options.plot_vol;

% Simulate trajectories and plot them on top
if plot_repr
    opt_sim = [];
    opt_sim.dt = 0.01;
    opt_sim.i_max = 5000;
    opt_sim.tol = 0.005;
    opt_sim.plot = 0;
    [x_sim, ~] = Simulation(x0_all ,[],ds_fun, opt_sim); 
else
	x_sim = [];
	hr = [];
end

% Plot Streamlines
if M == 2
    [hd] = scatter(Xi_ref(1,:),Xi_ref(2,:),10,[1 0 0],'filled'); hold on
    limits = axis;
    limits_ = limits + [-0.015 0.015 -0.015 0.015];
    [hs] = plot_ds_model(fig1, ds_fun, [0 0]', limits_,'medium'); hold on;
    axis(limits_)
    box on
    grid on
    xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
    ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);
    
    % Plot simulated trajectories
    if plot_repr
        [hr] = scatter(x_sim(1,:),x_sim(2,:),10,[0 0 0],'filled'); hold on
    end
    
elseif M == 3
    
    % Plot Demonstrations in red
    [hd] = plot3(Xi_ref(1,:),Xi_ref(2,:),Xi_ref(3,:),'r.','markersize',10); hold on;
       
    
    % Compute Start Locations for Streamlines     
    start_pnts =  sample_initial_points(x0_all, nb_pnts, init_type, plot_volume);
    limits = axis;
    
    % Plot Streamlines in blue
    [hs] = plot_ds_model_3D(fig1, ds_fun, [0;0;0], limits, start_pnts, 'low'); hold on;
    
    % Simulate trajectories and plot them on top
    if plot_repr
        [hr] = plot3(x_sim(1,:),x_sim(2,:),x_sim(3,:),'k.','markersize',10); hold on;
    end
    axis equal
end

end