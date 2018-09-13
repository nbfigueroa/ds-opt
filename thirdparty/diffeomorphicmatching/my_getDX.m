function [v, phiX] = my_getDX(A, ftauJ, X, phiZero, stepSize, epsilonTarg_2)

% Compute inverse Jacobian and Phi_inv_x
dim = size(X,1);
[phiX, Ji] = ftauJ(X);
lX = length(phiX);

%Relative dist for the Dynamic
phiXR = phiX - repmat(phiZero, 1, lX);

v = A*phiXR;
xn = sum(phiXR.^2,1);
vn = sqrt(sum(v.^2,1));
ind = xn > epsilonTarg_2;

% Hack to make the velocities at the attractor 0
if max(ind)
    v(:,ind) = v(:,ind).*repmat(stepSize./vn(ind),dim,1);
end
if ~min(ind)
    v(:,~ind) = v(:,~ind).*repmat(stepSize./vn(~ind),dim,1).*repmat(xn(~ind)./epsilonTarg_2,dim,1);
end

% Multiply jacobian times the estimated velocitoes
v = reshape(mtimesx(Ji, reshape(v,dim,1,lX)),dim,lX);

end

