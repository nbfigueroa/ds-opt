function [Xinit, target, targetV, source, sourceV, allTarget, allTargetV, allSource, allSourceV] = normalizeTrajectories(Xinit, target, targetV, source, sourceV, allTarget, allTargetV, allSource, allSourceV)
[dim, lX] = size(target);
varTarg = zeros(dim,1);


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