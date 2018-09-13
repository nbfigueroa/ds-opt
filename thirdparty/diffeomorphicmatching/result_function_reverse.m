function pt = result_function_reverse(centers, targets, coefs, division_coef, nb_iteration, pt, conv_crit)
%% Apply the inverse transformation of the transformation given by [centers, targets, coefs, division_coef, nb_iteration]
% to the points in pt
% The inverse transformation has no closed form, so we apply a newton-euler
% based minimization until conv_crit is reached

if nargin <= 6 || isempty(conv_crit)
    conv_crit = 1e-12;
end

if numel(division_coef) == 1
    division_coef = repmat(division_coef,1,nb_iteration);
end

[dim,lX] = size(pt);

for j = nb_iteration:-1:1      
    V = repmat(-(targets(:,j) - centers(:,j))/division_coef(j), 1, lX);
    b = pt - repmat(centers(:,j), 1, lX);
    
    lambda = zeros(1,lX);
    
    r = radial(centers(:,j), pt + repmat(lambda,dim,1) .* V, coefs(j));
    
    g = sum(2*(b+repmat(lambda,dim,1) .* V).*V, 1);
    deriv = -1 - coefs(j).^2 .* g .* r;
    errValue = r - lambda;    
    ind = abs(errValue) > conv_crit;

    while sum(ind)
        lambda(ind) = max(lambda(ind) - errValue(ind)./deriv(ind), 0);
        lambda(ind) = min(lambda(ind), 1);
        r(ind) = radial(centers(:,j), pt(:,ind) + repmat(lambda(ind),dim,1) .* V(:,ind), coefs(j));
        g(ind) = sum(2*(b(:,ind)+repmat(lambda(ind),dim,1) .* V(:,ind)).*V(:,ind), 1);
        deriv(ind) = -1 - coefs(j).^2 .* g(ind) .* r(ind);
        errValue(ind) = r(ind) - lambda(ind);
        ind = abs(errValue) > conv_crit;
    end    
    pt  = pt + repmat(lambda,dim,1) .* V;
end

end

function result = radial(center, pt, coef) 
result = exp(-coef.^2 .* sum( (pt - repmat(center,1,size(pt,2))).^2, 1) );
end