function h = plotForwardLyapFun2D( nFig, fun_reverse, xLim, N, ccolour, numPoints, LineWidth )
%PLOTLYAPFUN Summary of this function goes here
%   Detailed explanation goes here

%Options
if iscell(nFig)
    figure(nFig{1}); hold all;
    %set(gcf(), 'CurrentHandle', nFig{2})
    axes(nFig{2});
else
    figure(nFig); hold all;
end
axisOld = axis();

if nargin <= 2 || isempty(xLim)
    xLim = axis();
end

if nargin <= 3 || isempty(N)
    N=20;
end

if nargin <= 4 || isempty(ccolour)
    ccolor = jet(N);
end

if nargin <= 5 || isempty(numPoints)
    numPoints = 600;
end

if nargin <= 6 || isempty(LineWidth)
    LineWidth = 1;
end

Xmax = sqrt(max(xLim(1:2).^2)+max(xLim(3:4).^2));
%Try to guess a reasonable spacing
radius = linspace(sqrt(Xmax/N), sqrt(Xmax*(N-1)/N), N).^2;
numPoints = min(ceil(numPoints + numPoints/2.*radius.^2), 5000);
indPoints = [0,cumsum(numPoints)];

allPoints = zeros(2,indPoints(end));

for k = 1:N
    allV = repmat([1;0],1,numPoints(k));
    allAlpha = linspace(0,2*pi, numPoints(k));
    allAlphaS = sin(allAlpha);
    allAlphaC = cos(allAlpha);
    allR = zeros(2,2*numPoints(k));
    allR(1,1:2:end-1) = allAlphaC;
    allR(1,2:2:end) = -allAlphaS;
    allR(2,1:2:end-1) = allAlphaS;
    allR(2,2:2:end) = allAlphaC;
    allR = reshape(allR,2,2,numPoints(k));
    allPoints(:, 1+indPoints(k):indPoints(k+1)) = reshape(mtimesx(allR, reshape(allV.*radius(k), 2, 1, numPoints(k))), 2, numPoints(k));
end

allPoints = fun_reverse(allPoints);
h = cell(1,N);
for k = 1:N
    h{k} = plot(allPoints(1, 1+indPoints(k):indPoints(k+1)), allPoints(2, 1+indPoints(k):indPoints(k+1)), 'color', ccolor(k,:), 'LineWidth', LineWidth); 
end

axis(axisOld);

end
