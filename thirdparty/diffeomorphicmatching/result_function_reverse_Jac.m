function [pt, JacInv, Jac] = result_function_reverse_Jac(centers, targets, coefs, division_coef, nb_iteration, pt, JacInv, conv_crit)
%% Apply the inverse transformation of the transformation given by [centers, targets, coefs, division_coef, nb_iteration]
% to the points in pt
% The inverse transformation has no closed form, so we apply a newton-euler
% based minimization until conv_crit is reached

%Also returns the Jacobian of the transformation in each point

%% ATTENTION due to how the jacobian is calculated, by default to 'forward' jacobian [JacInv] is constructed. The jacobian of the reverse function [Jac] is only given on demand

if nargin <= 7 || isempty(conv_crit)
    conv_crit = 1e-12;
end

if numel(division_coef) == 1
    division_coef = repmat(division_coef,1,nb_iteration);
end

[dim,lX] = size(pt);

if nargin <= 6 || isempty(Jac)
    JacInv = repmat(eye(dim),1,1,lX);
end

Id = repmat(eye(dim),1,1,lX);

coefsS = coefs.^2;

for j = nb_iteration:-1:1      
    V = repmat(-(targets(:,j) - centers(:,j))/division_coef(j), 1, lX);
    thisCenter = repmat(centers(:,j), 1, lX);
    b = pt - thisCenter;
    
    lambda = zeros(1,lX);
    
    r = exp(-coefsS(j) .* sum( (pt + repmat(lambda,dim,1) .* V - thisCenter).^2, 1) );

    g = sum(2*(b+repmat(lambda,dim,1) .* V).*V, 1);
    deriv = -1 - coefsS(j) .* g .* r;
    errValue = r - lambda;    
    ind = abs(errValue) > conv_crit;

    while sum(ind)
        lambda(ind) = max(lambda(ind) - errValue(ind)./deriv(ind), 0);
        lambda(ind) = min(lambda(ind), 1);
        r(ind) = exp(-coefsS(j) .* sum( ((pt(:,ind) + repmat(lambda(ind),dim,1) .* V(:,ind)) - thisCenter(:,ind)).^2, 1) );
        g(ind) = sum(2*(b(:,ind)+repmat(lambda(ind),dim,1) .* V(:,ind)).*V(:,ind), 1);
        deriv(ind) = -1 - coefsS(j) .* g(ind) .* r(ind);
        errValue(ind) = r(ind) - lambda(ind);
        ind = abs(errValue) > conv_crit;
    end    
    pt  = pt + repmat(lambda,dim,1) .* V;
    
    %Attention to in which order the jacobians are multiplied
    thisDX = pt - thisCenter;
    %mtimesx has some overhead, so use regular matrix multiplication if
    %there is only one point in pt
    if lX == 1
        JacInv = JacInv * (Id+V*(thisDX'*lambda.*(2*coefsS(j))));
    else
        JacInv = mtimesx( JacInv, Id-bsxfun(@times, mtimesx( reshape(-1*V,dim,1,lX), reshape(thisDX,1,dim,lX), 'speed' ), reshape(lambda.*(2*coefsS(j)),1,1,lX)), 'speed' );
    end
end

if nargout > 2
    Jac = JacInv;
    for k = 1:lX
       Jac(:,:,k) = inv(JacInv(:,:,k)); 
    end
else
    Jac = [];
end

end

% function result = radial(center, pt, coef) 
% result = exp(-coef.^2 .* sum( (pt - repmat(center,1,size(pt,2))).^2, 1) );
% end