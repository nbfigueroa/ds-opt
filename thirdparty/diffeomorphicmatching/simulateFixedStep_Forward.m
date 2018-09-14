function [ PhiXsol, Xsol ] = simulateFixedStep_Forward( PhiX, ftauJ, X0, varargin )
%% Simulate the resulting movement

[dim,lX] = size(X0);
epsilonOrig = 0.05;
epsilonPhi = 0.01;
maxStep = 1001;
stepSizeFac = 1;
offsetEnd = zeros(dim,1);

EIG0 = -diag([1,2.*ones(1,dim-1)]);


% Converging linear Compliant DS
% A_c tracking Linear DS 
Y0 = PhiX(:,1);
Y0 = Y0./norm(Y0);
XZ0 = null(Y0).';
R0 = [Y0;XZ0].';
A = R0*EIG0*R0.';


epsilonTarg = epsilonOrig*norm(Y0);
epsilonTarg_2 = epsilonTarg^2;
stepSize = mean(sqrt(sum( diff(PhiX,[],2).^2,1)))*(size(PhiX,2))/maxStep*stepSizeFac;
epsilonTarg1 = epsilonPhi*sqrt(mean(sum(X0.^2,1)));
epsilonTarg1_2 = epsilonTarg1^2;



% Tracking linear Compliant DS
% A_t tracking Linear DS 
% y1 = 1;
% y2 = -Y0(1)/Y0(2);
% y = [y1;y2];
% Q = [y./norm(y),Y0./norm(Y0)];
% L = [-20 0 ; 0 -1];
% A = Q*L*Q';


phiZero = PhiX(:,end);
fDyn = @(X) getDX(A, ftauJ, X, phiZero, stepSize, epsilonTarg_2, dim, lX, 0);

PhiXsol = zeros(dim,lX,maxStep);
Xsol = zeros(dim,lX,maxStep);
Xsol(:,:,1) = X0;

for k=1:maxStep-1
    
    %Check break condition
    if min(  sum(Xsol(:,:,k).^2,1) < epsilonTarg1_2 )
        PhiXsol = PhiXsol(:,:,1:k);
        Xsol = Xsol(:,:,1:k);
        break;
    end
    [dX, PhiXsol(:,:,k)] = fDyn(Xsol(:,:,k));
    Xsol(:,:,k+1) = Xsol(:,:,k)+dX;
end

PhiXsol = PhiXsol(:,:,1:end-1);
Xsol = Xsol(:,:,1:end-1);
Xsol = Xsol -  repmat(offsetEnd, 1, size(Xsol,2), size(Xsol,3));
end
