function [ PhiXsol, Xsol ] = simulateFixedStep_Forward_Ideal( PhiX, fun, fun_reverse, X0, varargin )
%% Simulate the resulting trajectory in an idealized way:
% Get the straight line from the initial point to the origin and then
% transform it using the forward transformation.
% Note that this results in the exact same orbit as integrating the
% non-linear dynamics in the "demonstration" space in the continuous
% time case

[dim,lX] = size(X0);

EIG0 = -diag([1,2.*ones(1,dim-1)]);

maxStep = 1001;
stepSizeFac = 1;
offsetEnd = zeros(dim,1);

for k = 1:2:length(varargin)
    eval([varargin{k},'=', varargin{k+1}]);
end

Y0 = PhiX(:,1);
stepSize = mean(sqrt(sum( diff(PhiX,[],2).^2,1)))*(size(PhiX,2))/maxStep*stepSizeFac;

% Converging linear Compliant DS
% A_c tracking Linear DS 
Y0 = Y0./norm(Y0);
XZ0 = null(Y0).';
R0 = [Y0;XZ0].';
A = R0*EIG0*R0.';

phiZero = PhiX(:,end);

%'Ideal' integration ( so X is integrated in the transformed space and then
%the total is transformed

PhiXsol = zeros(dim, lX, maxStep-1);
%Relative distance for dynamic
PhiXsol(:,:,1) = fun_reverse(X0) - repmat(phiZero,1,lX);

for k = 1:maxStep-2
    v = A*PhiXsol(:,:,k);
    vn = sqrt(sum(v.^2,1));
    ind = vn > 3*stepSize;
    if max(ind)
        v(:,ind) = v(:,ind).*repmat(stepSize./vn(ind),dim,1);
    end
    if ~min(ind)
        v(:,~ind) = v(:,~ind).*0.2;
    end
    PhiXsol(:,:,k+1) = PhiXsol(:,:,k)+v;
end

%Total distance
PhiXsol = PhiXsol + repmat(phiZero,1,lX,(maxStep-1));
Xsol = reshape(fun(reshape(PhiXsol,dim,lX*(maxStep-1))), dim,lX,maxStep-1) - repmat(offsetEnd, 1, lX, maxStep-1);
end


