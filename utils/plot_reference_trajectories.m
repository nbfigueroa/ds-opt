function [h_data, h_att, h_vel] = plot_reference_trajectories(Data, att_g, att_l, vel_sample)

h_data = plot(Data(1,:),Data(2,:),'r.','markersize',10); hold on;
h_att = [];
% h_att = [h_att, scatter(att_g(1),att_g(2),150,[0 0 0],'d','Linewidth',2)]; hold on;
% 
% for k=1:size(att_l,2)    
%     h_att = [h_att, scatter(att_l(1,k),att_l(2,k),150,[0 0 0],'+','Linewidth',2)]; hold on;
% end
% Plot Velocities of Reference Trajectories
vel_points = Data(:,1:vel_sample:end);
U = zeros(size(vel_points,2),1);
V = zeros(size(vel_points,2),1);
for i = 1:size(vel_points, 2)
    dir_    = vel_points(3:end,i)/norm(vel_points(3:end,i));
    U(i,1)   = dir_(1);
    V(i,1)   = dir_(2);
end
h_vel = quiver(vel_points(1,:)',vel_points(2,:)', U, V, 0.5, 'Color', 'k', 'LineWidth',2); hold on;

end