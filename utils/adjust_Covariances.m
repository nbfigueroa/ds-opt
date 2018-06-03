function [Sigma] = adjust_Covariances(Sigma, tot_scale_fact, rel_scale_fact)


% Check Relative Covariance Matrix Eigenvalues
est_K = size(Sigma,3);
Vs = zeros(2,2,est_K);Ls = zeros(2,2,est_K);
p1_eig = []; p2_eig = [];
for k=1:est_K                    
    [Vs(:,:,k), Ls(:,:,k)] = eig(Sigma(:,:,k));           
    if ~issorted(diag(Ls(:,:,k)))
        [L_,ids] = sort(diag(Ls(:,:,k)))
        Ls(:,:,k) = diag(L_);
        Vs(:,:,k) = Vs(:,ids,k);
    end
    Ls(:,:,k) = tot_scale_fact*Ls(:,:,k);
    lambda_1 = Ls(1,1,k);
    lambda_2 = Ls(2,2,k);
    p1_eig = [p1_eig lambda_1];
    p2_eig = [p2_eig lambda_2];
    Sigma(:,:,k)    = Vs(:,:,k) * Ls(:,:,k) * Vs(:,:,k)';
end

% Scale Sigma's to increase the influence of the local linear dynamics
cov_ratios = p1_eig./p2_eig;
p2_scaled = p2_eig./max(p2_eig);
for k=1:est_K
%     if (cov_ratios(k) < scale_fact && p2_scaled(k) < 0.75)
      if (cov_ratios(k) < rel_scale_fact)
        lambda_1 = p1_eig(k); lambda_2 = p2_eig(k);
        lambda_1_ = lambda_1 + lambda_2*(rel_scale_fact-cov_ratios(k));
        new_ratio = lambda_1_/lambda_2;
        Sigma(:,:,k)    = Vs(:,:,k) * diag([ lambda_1_ ; lambda_2 ]) * Vs(:,:,k)';
    end
    
end

% TODO: Here I have to do a re-estimation of the priors

% Re-group Sigma
ds_gmm.Sigma = Sigma;
    
end