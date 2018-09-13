function result = iterative_function(pt1, pt2, pt, division_coef, coef)
%Function applying the transformation corresponding to one kernel
%The translation vector V = (pt2-pt1)/division_coef
[dim,l] = size(pt);
result = pt + repmat((pt2 - pt1)/division_coef,1,l).*repmat(exp(-((coef^2) .* sum((pt - repmat(pt1,1,l)).^2,1))),dim,1);
end