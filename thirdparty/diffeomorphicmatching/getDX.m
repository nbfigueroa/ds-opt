function [v, phiX] = getDX(A, ftauJ, X, phiZero, stepSize, epsilonTarg_2, dim, lX, doInvJac)

[phiX, Ji] = ftauJ(X);
lX = length(phiX);

%Relative dist for the Dynamic
phiXR = phiX - repmat(phiZero, 1, lX);

v = A*phiXR;
xn = sum(phiXR.^2,1);
vn = sqrt(sum(v.^2,1));
ind = xn > epsilonTarg_2;

if max(ind)
    v(:,ind) = v(:,ind).*repmat(stepSize./vn(ind),dim,1);
end
if ~min(ind)
    v(:,~ind) = v(:,~ind).*repmat(stepSize./vn(~ind),dim,1).*repmat(xn(~ind)./epsilonTarg_2,dim,1);
end
if doInvJac
    for k = 1:lX
        v(:,k) = Ji(:,:,k)\v(:,k);
    end
else
    if lX > 1
        v = reshape(mtimesx(Ji, reshape(v,dim,1,lX)),dim,lX);
    else
        v = mtimesx(Ji, v);
    end
end

end

