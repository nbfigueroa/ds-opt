function [ PhiXsol, Xsol ] = my_simulateFixedStep_Forward( PhiX, ftauJ, X0, A )
%% Simulate the resulting movement

[dim,lX] = size(X0);
epsilonOrig = 0.05;
epsilonPhi = 0.01;
maxStep = 1001;
stepSizeFac = 1;
offsetEnd = zeros(dim,1);

epsilonTarg = epsilonOrig*norm(Y0);
epsilonTarg_2 = epsilonTarg^2;
stepSize = mean(sqrt(sum( diff(PhiX,[],2).^2,1)))*(size(PhiX,2))/maxStep*stepSizeFac;
epsilonTarg1 = epsilonPhi*sqrt(mean(sum(X0.^2,1)));
epsilonTarg1_2 = epsilonTarg1^2;

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
