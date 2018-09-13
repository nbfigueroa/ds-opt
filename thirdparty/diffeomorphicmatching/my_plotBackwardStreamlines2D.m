function h = my_plotBackwardStreamlines2D( nFig, mean_target_traj, ftauJ, A )
%% Plot the streamlines resulting from the linear dynamics and the transformation

if iscell(nFig)
    handle = figure(nFig{1}); hold all;
    axes(nFig{2});
else
    handle = figure(nFig); hold all;
end

attractor = mean_target_traj(:,end);
diff_ds = @(x) nonlin_diff_ds(x, A, ftauJ, mean_target_traj, attractor);

% Create Grid and Evaluate Dynamics
xLim = axis();
h = plot_ds_model(handle, diff_ds, attractor, xLim);



end

