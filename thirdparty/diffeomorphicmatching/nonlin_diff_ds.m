function [x_dot] = nonlin_diff_ds(x, ftauJ, mean_target_traj)

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
phiXR = phiX;

% A of Linear DS 
A = -eye(2);

% Compute Dynamics with original linear model
x_dot = A*phiXR;

% Hack to make the velocities at the attractor 0
xn = sum(phiXR.^2,1);
vn = sqrt(sum(x_dot.^2,1));
ind = xn > epsilonTarg_2;
if max(ind)
    x_dot(:,ind) = x_dot(:,ind).*repmat(stepSize./vn(ind),dim,1);
end
if ~min(ind)
    x_dot(:,~ind) = x_dot(:,~ind).*repmat(stepSize./vn(~ind),dim,1).*repmat(xn(~ind)./epsilonTarg_2,dim,1);
end

% Multiply jacobian times the estimated velocities for non-linear DS
x_dot = reshape(mtimesx(Ji, reshape(x_dot,dim,1,lX)),dim,lX);

end

