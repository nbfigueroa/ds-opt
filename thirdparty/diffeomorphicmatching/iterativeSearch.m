function [ centers, targets, coefs, division_coef, nb_iteration ] = iterativeSearch( source, target, varargin )
%ITERATIVESEARCH Searches for the parameters of the transformation
%With explicit end matching

[dim, lX] = size(source);

%Default options
nb_iteration = 150;
maxCoef = 50;
conv_crit = 1e-6;
doPlot = 0;
regularise = 5e-4;
deltaPlot =1;

points = [0,1,2,3,4,5,6]./6;
division_coefList = [2,2,2,1.5,1,1,1];
safeCoeffList = [0.7, 0.7, 0.7 0.65, 0.6, 0.6, 0.6];

for k = 1:2:length(varargin)
    eval([varargin{k},'=', varargin{k+1}, ';']);
end

%Initialize output
centers = zeros(dim,nb_iteration-1);
targets = zeros(dim,nb_iteration-1);
coefs = zeros(1,nb_iteration-1);

%Interpolate parameters for each step linearily
points = floor(points.*(nb_iteration-1));
[points, ii] = unique(points);
division_coefList = division_coefList(ii);
safeCoeffList = safeCoeffList(ii);
division_coef = interp1(points, division_coefList, 1:(nb_iteration-1), 'linear');
safeCoeff = interp1(points, safeCoeffList, 1:(nb_iteration-1), 'linear');

%% Start the search
R = source;
for j = 1:nb_iteration-1
    [themax, idmax] = max(sum((R - target).^2,1));    
    if themax < conv_crit
        nb_iteration = j-1;
        coefs = coefs(1:nb_iteration);
        centers = centers(:,1:nb_iteration);
        targets = targets(:,1:nb_iteration);
        division_coef = division_coef(:,1:nb_iteration);
        break;
    end
    
    if doPlot && mod(j, deltaPlot) == 0
        figure(); hold all;
        plot(R(1,:), R(2,:), '-b', 'linewidth', 2);
        v = (target(:,idmax)-R(:,idmax))./division_coef(j);
        plot([R(1,idmax), target(1,idmax)], [R(2,idmax), target(2,idmax)], '.-b');
        plot([R(1,idmax), R(1,idmax)+v(1)], [R(2,idmax), R(2,idmax)+v(2)], '.-k');
        plot(target(1,:), target(2,:), 'g');    
    end

    % Maximal coefficient allowed to ensure diffeomorphic properties
    coefmax = safeCoeff(j) * division_coef(j) * 1/(sqrt(2) * exp(-1/2) * norm(target(:,idmax) - R(:,idmax)));
    
    % Cost function is a mixture of the distance between the transformed
    % source points to target and a regularization value preventing
    % overdeformation
    f = @(x) norm(iterative_function(R(:,idmax), target(:,idmax), R, division_coef(j), x)-target)./lX + regularise*(x./coefmax)^2;
    
    %Single parameter bounded search
    [coef, fval] = fminbnd(f,0,min(coefmax, maxCoef));
    
    centers(:,j) = R(:,idmax);
    targets(:,j) = target(:,idmax);
    coefs(j) = coef;

    if doPlot && mod(j, deltaPlot) == 0
        lx = get(gca(), 'xlim');
        ly = get(gca(), 'ylim');
        PP = getPoints(coef);
        plot(PP(1,:)+R(1,idmax), PP(2,:)+R(2,idmax), 'g');
        xlim(lx);
        ylim(ly);
        title(num2str(sum((sum(((iterative_function(R(:,idmax), target(:,idmax), R, division_coef(j), coef) - target)).^2,1)).^0.5)));
    end
    %Apply the iteration found in this step
    R = iterative_function(R(:,idmax), target(:,idmax), R, division_coef(j), coef);
    if doPlot && mod(j, deltaPlot) == 0
        hh = plot(R(1,:), R(2,:), '--');
        plot(R(1,:), R(2,:), '-r', 'linewidth', 2);
    end
        
end

%Do the final matching step such that the transformation at the origin is
%the Ientity

themax = sum((R(:,end) - target(:,end)).^2,1);
idmax = size(R, 2);
division_coefEnd = 1;

if norm(themax) >= 1e-6
    f = @(x) norm(iterative_function(R(:,idmax), target(:,idmax), R, division_coef(j), x)-target)./lX + regularise*(x./coefmax)^2;
    coefmax = safeCoeff(end) * division_coefEnd * 1/(sqrt(2) * exp(-1/2) * norm(target(:,idmax) - R(:,idmax)));
    coef = fminbnd(f,0,min([coefmax, maxCoef, 1e6]));
    centers(:,end+1) = R(:,idmax);
    targets(:,end+1) = target(:,idmax);
    coefs(end+1) = coef;
    division_coef(end+1) = division_coefEnd;
    safeCoeff(end+1) = safeCoeff(end);
    R = iterative_function(R(:,idmax), target(:,idmax), R, division_coefEnd, coef);
end
nb_iteration = length(division_coef);

end

function p = getPoints(coef)
%Helper fonction to plot region of influence
ang = [linspace(0, 2*pi*19/20, 20),0];
p=[cos(ang); sin(ang)].*(-1/coef^2*log(0.5));
end
