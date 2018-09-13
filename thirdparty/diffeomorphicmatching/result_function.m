function result = result_function(centers, targets, coefs, division_coef, nb_iteration, pt)
%% Applies the forward transformation given as [centers, targets, coefs, division_coef, nb_iteration]
% to the points in pt

[dim,l] = size(pt);
if numel(division_coef) == 1
    division_coef = repmat(division_coef,1,nb_iteration);
end
for j = 1:nb_iteration
    pt = pt + repmat((targets(:,j) - centers(:,j))/division_coef(j),1,l).*repmat(exp(-(coefs(j)^2 .* sum((pt - repmat(centers(:,j),1,l)).^2,1))),dim,1);
end
result = pt;
end