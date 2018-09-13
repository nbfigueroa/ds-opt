function [Xinit, target, allTarget, allTargetV, allSource, allSourceV, nDemos, indD] = processDrawnDataset(data, use_demos, step_demos,dt)
% Process Drawn Dataset Function
dim =2;


% Pre-process demos so they have the same number of points
lengths = [];
for d=1:length(data)
    lengths = [lengths length(data{d})]; 
end
max_length = max(lengths);
for d=1:length(data)
    s = length(data{d})/max_length;
    data{d} = [interp1(1:length(data{d}),data{d}(1,:),1:s:length(data{d}),'spline');
                interp1(1:length(data{d}),data{d}(2,:),1:s:length(data{d}),'spline');
                interp1(1:length(data{d}),data{d}(3,:),1:s:length(data{d}),'spline');
                interp1(1:length(data{d}),data{d}(4,:),1:s:length(data{d}),'spline')];
    if length(data{d})== max_length
        data{d} = data{d}(:,1:end-1);
    end
end

allTarget = [];
allTargetV = [];
allSource = [];
allSourceV = [];
Xinit = [];
[~, lX] = size(data{1});
indD = 1:step_demos:lX;
if indD(end) ~= lX
    indD = [indD, lX];%#k
end


target = 0;
nDemos = length(use_demos);
for k = 1:nDemos
    demo = data{use_demos(k)};
    % Initial points of demonstrations
    Xinit = [Xinit, demo(1:2,1)];
    % Target trajectories
    target = target+demo(1:2,:);
    % All Target trajectories
    allTarget = [allTarget,   demo(1:2,indD)];
    allTargetV = [allTargetV, demo(3:4,indD)];
    
    % This source trajectories
    thisSource  = [linspace(demo(1,1), demo(1,end), lX);linspace(demo(2,1), demo(2,end), lX)];
    thisSourceV = [diff(thisSource./dt, [], 2), zeros(dim,1)];
    
    % All source trajectories
    allSource = [allSource, thisSource(:,indD)];%#ok
    allSourceV = [allSourceV, thisSourceV(:,indD)];%#ok
end
end