clear all; close all;

%All Lasa-shapes except multidemonstration examples
dataset_dir = '../../LASADataset/Dataset';
allShapes = dir(dataset_dir);
allShapes(1:2) = [];
figName = 'recorded_motions/';
suffix = '';
allOptions = {};
mkdir(figName);

%% Options
%Options concerning th simulation
%Set eigenvalues here
%'EIG0' defines the linear system. EIG0(2,2) specifies the convergence rate
%in the orthogonal direction
maxStep = 1001;
optsSim = {'maxStep', num2str(maxStep), 'stepSizeFac',' 1', 'EIG0', '-diag([1,1.5])','offsetEnd', '[0;0]' }; 

%Other options
allInOnePlot = false;
plotGridNotTransform = false;
step_demos = 1; %Step the points
use_demos = []; %Use only the demos with the appearing here; Use all if empty
doNormalize = 1; %Normalize the demos using the variance in each direction; This works better in general

figure(1); hold all;%general plot
figure(2); hold all;%Streamlines
figure(3); hold all;%Lyap function level-sets

%Options concerning the transformation
if doNormalize
    optsSearch = {'maxCoef', '10', 'nb_iteration','150', 'regularise', '5e-4', 'conv_crit', '1e-6', 'division_coefList', '[3,3,2.5,2.0,1.5,1.1,1.1]', 'safeCoeffList', '0.6*[1, 1, 1, 1, 1, 1, 1]', 'doPlot', '0'};
else
    optsSearch = {'maxCoef', '5', 'nb_iteration','150', 'regularise', '1e-3', 'conv_crit', '1e-6', 'division_coefList', '[3,3,2.5,2.0,1.5,1.1,1.1]', 'safeCoeffList', '0.55*[1, 1, 1, 1, 1, 1, 1]', 'doPlot', '0'};
end

% for f = 1:length(allShapes);
for f = 1:1;    

    %% Load data
    thisFile = allShapes(f).name;
    load([dataset_dir, thisFile]);
    if isempty(use_demos)
        use_demos = 1:length(demos);
    end
    fprintf('*** Loaded demonstrations from %s ***\n', allShapes(f).name);
    %% Prepare data
    target = 0;
    allTarget = [];
    allTargetV = [];
    allSource = [];
    allSourceV = [];
    Xinit = [];
    [dim, lX] = size(demos{1}.pos);
    indD = 1:step_demos:lX;
    if indD(end) ~= lX
        indD = [indD, lX];%#k
    end
    
    fprintf('- Preparing Data for Diffeomorphic Matching...');
    nDemos = 0;
    for k = 1:length(use_demos)
        nDemos = nDemos+1;
        demos{use_demos(k)}.pos(:,end);
        Xinit = [Xinit, demos{use_demos(k)}.pos(:,1)];%#ok
        target = target+demos{use_demos(k)}.pos;
        allTarget = [allTarget, demos{use_demos(k)}.pos(:,indD)];%#ok
        thisV = [diff(demos{k}.pos,[],2)./dt, zeros(dim,1)];
        allTargetV = [allTargetV, thisV(:,indD)];%#ok
        thisSource = [linspace(demos{use_demos(k)}.pos(1,1), demos{use_demos(k)}.pos(1,end), lX); linspace(demos{use_demos(k)}.pos(2,1), demos{use_demos(k)}.pos(2,end), lX)];
        thisSourceV = [diff(thisSource./dt, [], 2), zeros(dim,1)];
        allSource = [allSource, thisSource(:,indD)];%#ok
        allSourceV = [allSourceV, thisSourceV(:,indD)];%#ok
    end

    fprintf('done\n');
    % Use the mean of all demonstrations as target (The transformation maps
    % the points given in "source" onto the points in "target"
    target = target./nDemos;
    targetV = [diff(target,[],2),zeros(dim,1)]./dt;
    
    %Define the source: A straight line betwenn the initial and final point
    source = [linspace(target(1,1), target(1,end), lX); linspace(target(2,1), target(2,end), lX)];
    sourceV = [diff(source, [], 2).*(step_demos/dt), zeros(dim,1)];
    source = source(:,indD);
    sourceV = sourceV(:,indD);
    
    target = target(:,indD);
    targetV = targetV(:,indD);
    
    [dim, lX] = size(target);
    
    varTarg = zeros(dim,1);
    
    if doNormalize
        for k = 1:dim
            varTarg(k) = sqrt(var(allTarget(k,:)));
            target(k,:) = target(k,:)./varTarg(k);
            allTarget(k,:) = allTarget(k,:)./varTarg(k);%#ok
            targetV(k,:) = targetV(k,:)./varTarg(k);
            allTargetV(k,:) = allTargetV(k,:)./varTarg(k);%#ok
            source(k,:) = source(k,:)./varTarg(k);
            allSource(k,:) = allSource(k,:)./varTarg(k);%#ok
            sourceV(k,:) = sourceV(k,:)./varTarg(k);
            allSourceV(k,:) = allSourceV(k,:)./varTarg(k);%#ok
            Xinit(k,:) = Xinit(k,:)./varTarg(k);%#ok
        end
    end
    
    %Max for plotting
    XLimPlot = [min([1.2*allTarget, 0.8*allTarget, -.25*ones(dim,1)],[], 2), max([1.2*allTarget, 0.8*allTarget, .25*ones(dim,1)],[],2)];
    XLimPlot2 = [min([1.2*allSource, 0.8*allSource, -.25*ones(dim,1)],[], 2), max([1.2*allSource, 0.8*allSource, .25*ones(dim,1)],[],2)];
    
    %Search for the transformation parameters
    [ centers, targets, coefs, division_coef, nb_iteration ] = iterativeSearch( source, target, optsSearch{:} );

    %% Create all lambda function needed to evaluate the transformation
    % Forward transformation
    fun = @(pt) result_function(centers, targets, coefs, division_coef, nb_iteration, pt);
    % Backward transformation
    fun_reverse = @(pt) result_function_reverse(centers, targets, coefs, division_coef, nb_iteration, pt);
    % Backward transformation and calculation of the jacobian in point
    fun_reverse_jac = @(pt) result_function_reverse_Jac(centers, targets, coefs, division_coef, nb_iteration, pt);
 
    %General plot with the transformed sources and targets
    figure(1); hold all;
    if allInOnePlot
        subplot(5,6,f-2); hold all;
        axis tight
        axis off
    else
        clf;hold all;
    end

    plot(target(1,:),target(2,:), '-k', 'linewidth',2);
    invPhiTarg = fun_reverse(target);%Applying the inverse transformation to the mean target
    plot(invPhiTarg(1,:), invPhiTarg(2,:),'--b'); hold on;

    plot(source(1,:),source(2,:), '-b', 'linewidth',2);
    invPhiSrc = fun(source);
    plot(invPhiSrc(1,:), invPhiSrc(2,:),'r', 'linewidth',2); hold on;
    
    if plotGridNotTransform
        plotGrid(1, XLimPlot2, 10, 1000, fun);
    else
        phiAllSrc = fun(allSource);
        for k = 1:nDemos
            plot(allTarget(1,(k-1)*lX+1:k*lX),allTarget(2,(k-1)*lX+1:k*lX), '-c', 'linewidth',2);
            plot(phiAllSrc(1,(k-1)*lX+1:k*lX), phiAllSrc(2,(k-1)*lX+1:k*lX),'-m', 'linewidth',2); hold on;
        end
    end

    xlim(XLimPlot(1,:));
    ylim(XLimPlot(2,:));
      
    %Plot the resulting Lyapunov functions (transformation applied to
    %x²+y²==alpha²
    figure(2); hold all;
    if allInOnePlot
        subplot(5,6,f-2); hold all;
        axis tight
        axis off
    else
        clf; hold all;
    end
    
    xlim(XLimPlot(1,:));
    ylim(XLimPlot(2,:));
    hL = plotForwardLyapFun2D({gcf(), gca()}, fun, [], 35);
    for k = 1:nDemos
        plot(allTarget(1,(k-1)*lX+1:k*lX),allTarget(2,(k-1)*lX+1:k*lX), '-k');
    end

    %Plot the "dynamics"
    h = figure(3); hold all;
    if allInOnePlot
        subplot(5,6,f-2); hold all;
        axis tight
        axis off
    else
        clf;hold all;
    end
    xlim(XLimPlot(1,:));
    ylim(XLimPlot(2,:));
    hSL = plotBackwardStreamlines2D( {gcf(), gca()}, source, fun_reverse_jac, [], optsSim{:} );
    set(hSL,'LineWidth',1)
%     [phiX, X] = simulateFixedStep_Forward(source, fun_reverse_jac, Xinit, optsSim{:});
    [phiX, X] = simulateFixedStep_Forward_Ideal( source, fun, fun_reverse, Xinit, optsSim{:} );
    for k = 1:nDemos
        plot(allTarget(1,(k-1)*lX+1:k*lX),allTarget(2,(k-1)*lX+1:k*lX), '-k', 'LineWidth', 1);
    end
    for k = 1:size(X,2)
        plot(squeeze(X(1,k,:)), squeeze(X(2,k,:)), '-r', 'LineWidth', 1);
    end
    
    %Plot the velocity profiles
    h = figure(4); hold all;
    if allInOnePlot
        subplot(5,6,f-2); hold all;
        axis tight
        axis off
    else
        clf;hold all;
    end
    
    Tdemos = (0:(size(demos{1}.pos,2)-1)).*demos{1}.dt;
    dTsim = demos{1}.dt*size(demos{1}.pos,2)/maxStep;
    Tsim = (0:size(X,3)-2)*dTsim;
    for k = 1:nDemos
        if doNormalize
            plot(Tdemos, demos{k}.vel(1,:)./varTarg(1), '--k');
            plot(Tdemos, demos{k}.vel(2,:)./varTarg(2), '--r');
        else
            plot(Tdemos, demos{k}.vel(1,:), '--k');
            plot(Tdemos, demos{k}.vel(2,:), '--r');
        end
        thisVx = diff(squeeze(X(1,k,:)),1)./dTsim;
        thisVy = diff(squeeze(X(2,k,:)),1)./dTsim;
        plot(Tsim, thisVx, '-k');
        plot(Tsim, thisVy, '-r');
    end
    
    
    if ~allInOnePlot
        [a,b,c] = fileparts(thisFile);
%         saveas(figure(1), [figName, b,'_1', suffix,'.svg']);
        print(figure(1), '-dpng', '-r300', [figName, b, '_1', suffix,'.png']);
%         saveas(figure(2), [figName, b, '_2', suffix,'.svg']);
        print(figure(2), '-dpng', '-r300', [figName, b, '_2', suffix,'.png']);
%         saveas(figure(3), [figName, b, '_3', suffix,'.svg']);
        print(figure(3), '-dpng', '-r300', [figName, b, '_3', suffix,'.png']);
%         saveas(figure(4), [figName, b, '_4', suffix,'.svg']);
        print(figure(4), '-dpng', '-r300', [figName, b, '_4', suffix,'.png']);
    end

drawnow
end
if allInOnePlot
    b = 'allPlots';
    saveas(figure(1), [figName, b,'_1', suffix,'.svg']);    
    print(figure(1), '-dpng', '-r300', [figName, b, '_1', suffix,'.png']);
    saveas(figure(2), [figName, b, '_2', suffix,'.svg']);
    print(figure(2), '-dpng', '-r300', [figName, b, '_2', suffix,'.png']);
    saveas(figure(3), [figName, b, '_3', suffix,'.svg']);
    print(figure(3), '-dpng', '-r300', [figName, b, '_3', suffix,'.png']);
end

