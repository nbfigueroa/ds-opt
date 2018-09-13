function [Xinit, target, allTarget, allTargetV, allSource, allSourceV, nDemos, indD] = processLASADataset(demos, use_demos, step_demos, dt)
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

target = 0;
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

end