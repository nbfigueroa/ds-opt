function [hd, hs, hr, x_sim] = visualizeEstimatedDS(Xi_ref, ds_fun, varargin)
fig1 = figure('Color',[1 1 1]);
[hd] = scatter(Xi_ref(1,:),Xi_ref(2,:),10,[1 0 0],'filled'); hold on
limits = axis;
limits_ = limits + [-0.015 0.015 -0.015 0.015];
[hs] = plot_ds_model(fig1, ds_fun, [0 0]', limits_,'medium'); hold on;
axis(limits_)
box on
grid on
xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);

if nargin > 2
    plot_repr = varargin{1};
    x0_all    = varargin{2};
    % Simulate trajectories and plot them on top
    if plot_repr
        opt_sim = [];
        opt_sim.dt = 0.01;
        opt_sim.i_max = 3000;
        opt_sim.tol = 0.1;
        opt_sim.plot = 0;
        [x_sim xd_sim]=Simulation(x0_all ,[],ds_fun, opt_sim);
        [hr] = scatter(x_sim(1,:),x_sim(2,:),10,[0 0 0],'filled'); hold on
    end
else
    x_sim = [];
end
end