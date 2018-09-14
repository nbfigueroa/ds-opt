function h = plotBackwardStreamlines2D( nFig, PhiX, ftauJ, xLim, varargin )
%% Plot the streamlines resulting from the linear dynamics and the transformation

if iscell(nFig)
    figure(nFig{1}); hold all;
    %set(gcf(), 'CurrentHandle', nFig{2})
    axes(nFig{2});
else
    figure(nFig); hold all;
end
axisOld = axis();

if nargin <= 3 || isempty(xLim)
    xLim = axis();
end

dim = size(PhiX,1);

EIG0 = -diag([1,2.*ones(1,dim-1)]);
epsilonOrig = 0.05;
numPoints = 300;
doStreams = 1;

for k = 1:2:length(varargin)
    eval([varargin{k}, '=', varargin{k+1}]);
end

Y0 = PhiX(:,1);
epsilonTarg = epsilonOrig*norm(Y0);
epsilonTarg_2 = epsilonTarg^2;
maxStep = 1001;
stepSizeFac = 1;

stepSize = mean(sqrt(sum( diff(PhiX,[],2).^2,1)))*(size(PhiX,2))/maxStep*stepSizeFac;

Y0 = Y0./norm(Y0);
XZ0 = null(Y0).';
R0 = [Y0;XZ0].';
A = R0*EIG0*R0.';

[X, Y] = meshgrid(linspace(xLim(1), xLim(2), numPoints), linspace(xLim(3), xLim(4), numPoints));
allX = [X(:)'; Y(:)'];

phiZero = PhiX(:,end);

fDyn = @(X) getDX(A, ftauJ, X, phiZero, stepSize, epsilonTarg_2, dim, size(allX,2), 0);

allV = fDyn(allX);

if doStreams
    h = streamslice(reshape(allX(1,:), numPoints, numPoints), reshape(allX(2,:), numPoints, numPoints), reshape(allV(1,:), numPoints, numPoints), reshape(allV(2,:), numPoints, numPoints));
else
    h = quiver(allX(1,:), allX(2,:), allV(1,:), allV(2,:)); %#ok
end

axis(axisOld);

end

