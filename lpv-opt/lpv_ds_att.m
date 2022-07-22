function [x_dot] = lpv_ds_att(x, att_orig, att_new, Rot_att, ds_gmm, A_k)

% Auxiliary Variables
[N,M] = size(x);
K = length(ds_gmm.Priors);

% Works only for linear motion changes in (x,y,z) for rotating the DS we
% need to compute an isomorphism (basically full H matrix between original
% and new attractors)
att_offset = -att_new + att_orig;

% Posterior Probabilities per local DS
beta_k_x = posterior_probs_gmm(x + att_offset,ds_gmm,'norm');

% Output Velocity
x_dot = zeros(N,M);
for i = 1:size(x,2)
    % Estimate Global Dynamics component as LPV
    if K > 1
        f_k = zeros(N,K);
        for k=1:K
            f_k(:,k) = beta_k_x(k,i) * A_k(:,:,k)*inv(Rot_att)*(x(:,i) + att_offset - att_orig);
        end
        f_k = sum(f_k,2);
    else
        % Estimate Global Dynamics component as Linear DS
        f_k = A_k*inv(Rot_att)*(x(:,i) + att_offset - att_orig);
    end
    x_dot(:,i) = Rot_att*f_k;          
end
end



