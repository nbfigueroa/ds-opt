function [xdot] = diffeomorphic_ds(x, EIG0, PhiX, ftauJ)

N = size(x,2);
M = size(PhiX,1);
epsilonOrig = 0.05;
Y0 = PhiX(:,1);
epsilonTarg = epsilonOrig*norm(Y0);
epsilonTarg_2 = epsilonTarg^2;
Y0 = Y0./norm(Y0);
XZ0 = null(Y0).';
R0 = [Y0;XZ0].';
A = R0*EIG0*R0.';

phiZero = PhiX(:,end);
fDyn = @(X) getDX(A, ftauJ, X, phiZero, 1, epsilonTarg_2, M, N, 0);
xdot = fDyn(x);

if N == 1
    xdot = xdot(:,1);
end

end