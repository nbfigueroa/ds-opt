function plotGrid(nFig, xLim, N, N2, fun, LineStyle)
%% Plot the image of a regular grid by fun (can be a transformation or its inverse)
if nargin <= 0 || isempty(nFig)
    nFig = gcf();
end
figure(nFig); hold all;
if nargin <= 1 || isempty(xLim)
    xLim = (reshape(axis(),2,2))';
end
if nargin <= 2 || isempty(N)
    N = 20;
end
if nargin <= 3 || isempty(N2)
    N2 = 3000;
end
if nargin <= 4 || isempty(fun)
    fun = @(x) x;
end
if nargin <= 5 || isempty(LineStyle)
    LineStyle = {'-k', '-k'};
end

points = [];
for k =1:N
    points = [points, [(xLim(1,1)+(xLim(1,2)-xLim(1,1))*(k-1)/(N-1))*ones(1,N2); linspace(xLim(2,1), xLim(2,2), N2)], ...
                      [linspace(xLim(1,1), xLim(1,2), N2); (xLim(2,1)+(xLim(2,2)-xLim(2,1))*(k-1)/(N-1))*ones(1,N2)]];
end
points = fun(points);
for k = 1:N
    plot(points(1,1+2*N2*(k-1):2*N2*(k-1)+N2), points(2,1+2*N2*(k-1):2*N2*(k-1)+N2), LineStyle{1},'linewidth', 0.25);
    plot(points(1,1+N2+2*N2*(k-1):N2+2*N2*(k-1)+N2), points(2,1+N2+2*N2*(k-1):N2+2*N2*(k-1)+N2), LineStyle{2}, 'linewidth', 0.25);
end

end