function[h_vel] = visualizeEstimatedVelocities(Data, ds_fun)

[M, N] = size(Data);
M = M/2;
h_vel = figure('Color',[1 1 1]);
xd_dot = [];

% Simulate velocities from same reference trajectory
for i=1:N
    xd_dot_ = ds_fun(Data(1:M,i));    
    % Record Trajectories
    xd_dot = [xd_dot xd_dot_];        
end

% Plot Demonstrated Velocities vs Generated Velocities
if M == 2
    plot(Data(3,:)', '.-','Color',[0 0 1], 'LineWidth',1); hold on;
    plot(Data(4,:)', '.-','Color',[1 0 0], 'LineWidth',1); hold on;
    plot(xd_dot(1,:)','--','Color',[0 0 1], 'LineWidth', 1); hold on;
    plot(xd_dot(2,:)','--','Color',[1 0 0], 'LineWidth', 1); hold on;
    legend({'$\dot{\xi}^{ref}_{1}$','$\dot{\xi}^{ref}_{2}$','$\dot{\xi}^{d}_{1}$','$\dot{\xi}^{d}_{2}$'}, 'Interpreter', 'LaTex', 'FontSize', 15)    
elseif M == 3
    plot(Data(4,:)', '.-','Color',[0 0 1], 'LineWidth',1); hold on;
    plot(Data(5,:)', '.-','Color',[1 0 0], 'LineWidth',1); hold on;
    plot(Data(6,:)', '.-','Color',[0 1 0], 'LineWidth',1); hold on;
    plot(xd_dot(1,:)','--','Color',[0 0 1], 'LineWidth', 1); hold on;
    plot(xd_dot(2,:)','--','Color',[1 0 0], 'LineWidth', 1); hold on;
    plot(xd_dot(3,:)','--','Color',[0 1 0], 'LineWidth', 1); hold on;
    legend({'$\dot{\xi}^{ref}_{1}$','$\dot{\xi}^{ref}_{2}$','$\dot{\xi}^{ref}_{3}$',...,
        '$\dot{\xi}^{d}_{1}$','$\dot{\xi}^{d}_{2}$', '$\dot{\xi}^{d}_{3}$'}, 'Interpreter', 'LaTex', 'FontSize', 15)            
end
grid on;
title('Real vs Estimated Velocities', 'Interpreter', 'LaTex', 'FontSize', 15)
end