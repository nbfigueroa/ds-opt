function [xdot] = nonlin_diff_ds(x, A, ftauJ, mean_target_traj, attractor)

% Parameters for DS
epsilonOrig = 0.05;
Y0 = mean_target_traj(:,1);
epsilonTarg = epsilonOrig*norm(Y0);
epsilonTarg_2 = epsilonTarg^2;
maxStep = 1001;
stepSize = mean(sqrt(sum( diff(mean_target_traj,[],2).^2,1)))*(size(mean_target_traj,2))/maxStep;

% Compute inverse Jacobian and Phi_inv_x
dim = size(x,1);
[phiX, Ji] = ftauJ(x);
lX = length(phiX);

%Relative dist to attractor for the Dynamic
phiXR = phiX - repmat(attractor, 1, lX);

% Compute Dynamics with original linear model
xdot = A*phiXR;

% Hack to make the velocities at the attractor 0
xn = sum(phiXR.^2,1);
vn = sqrt(sum(xdot.^2,1));
ind = xn > epsilonTarg_2;
if max(ind)
    xdot(:,ind) = xdot(:,ind).*repmat(stepSize./vn(ind),dim,1);
end
if ~min(ind)
    xdot(:,~ind) = xdot(:,~ind).*repmat(stepSize./vn(~ind),dim,1).*repmat(xn(~ind)./epsilonTarg_2,dim,1);
end

% Multiply jacobian times the estimated velocities for non-linear DS
xdot = reshape(mtimesx(Ji, reshape(xdot,dim,1,lX)),dim,lX);

end