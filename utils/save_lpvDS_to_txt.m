function save_lpvDS_to_txt(DS_name, pkg_dir,  ds_gmm, A_k, b_k, att)

model_dir = strcat(pkg_dir,'/models/',DS_name, '/');
mkdir(model_dir)

% GMM parameters
Priors = ds_gmm.Priors;
Mu     = ds_gmm.Mu;
Sigma  = ds_gmm.Sigma;

% Writing Priors
dlmwrite(strcat(model_dir,'Priors'), Priors,'newline','pc','Delimiter',' ','precision','%.6f');

% Writing Mu
dlmwrite(strcat(model_dir,'Mu'), Mu,'newline','pc','Delimiter',' ','precision','%.6f');

% Writing Sigma
for i=1:length(Priors)   
    dlmwrite(strcat(model_dir,'Sigma'), Sigma(:,:,i),'newline','pc','-append','Delimiter',' ','precision','%.6f');
end

% Writing A's
for i=1:length(Priors)   
    dlmwrite(strcat(model_dir,'A_k'), A_k(:,:,i),'newline','pc','-append','Delimiter',' ','precision','%.6f');
end

% Writing b_k's
dlmwrite(strcat(model_dir,'b_k'), b_k,'newline','pc','Delimiter',' ','precision','%.6f');

% Writing attractor
dlmwrite(strcat(model_dir,'Attractor'), att,'newline','pc','Delimiter',' ','precision','%.6f');

% Writing Dimensions
Dimensions = [length(Priors); size(Mu,1)];
dlmwrite(strcat(model_dir,'Dimenions'), Dimensions,'newline','pc','Delimiter',' ','precision','%.6f');
