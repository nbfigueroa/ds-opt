function structGMM = learn_coupling(Data, master_dim, slave_dim, coupling_func, K)
%
% This function learns the coupling model between the selected dimensions
% using the specified coupling function.
%
% Inputs -----------------------------------------------------------------
%
%
%   o Data:    A 2d x N_Total matrix containing all the d-dimensional 
%              demonstration data points.
%              Rows 1:d corresponds to trajectories and the rows d+1:2d
%              are their first time derivatives. Each column of Data stands
%              for a datapoint. All demonstrations are put next to each other 
%              along the second dimension. For example, if we have 3 demos
%              D1, D2, and D3, then the matrix Data is:
%                               Data = [[D1] [D2] [D3]]
%
%   o master_dim:      Indices of the dimensions to be selected as the master
%                      sub-system
%
%   o slave_dim:       Indices of the dimensions to be selected as the slave
%                      sub-system
%
%   o coupling_func:  Function handle representing a transformation of the
%                     master state which will be coupled with the slave
%
%   o K:              Number of gaussians to be used for modelling
%
% Outputs ----------------------------------------------------------------
%
%   o structGMM:      A struct containing the coupling model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             Copyright (c) 2013 A. Shukla, LASA, EPFL,               %%%
%%%          CH-1015 Lausanne, Switzerland, http://lasa.epfl.ch         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The program is free for non-commercial academic use. Please contact the
% author if you are interested in using the software for commercial purposes.
% The software must not be modified or distributed without prior permission
% of the authors. Please acknowledge the authors in any academic publications
% that have made use of this code or part of it. Please use this BibTex
% reference:
% 
% Shukla, A. and Billard, A. " Coupled dynamical system based arm-hand 
% grasping model for learning fast adaptation strategies", Robotics and
% Autonomous Systems 2012.
%
%
% Please send your feedbacks or questions to:
%                           ashwini.shukla_at_epfl.ch

if(isempty(coupling_func))
    coupling_func = @(x) norm(x);
else if ~isa(coupling_func, 'function_handle')
        structGMM=[];
        disp('Invalid coupling function!');
        return;
    end
end

% Collecting the master and slave dimensions separately
master_data=[];
for i=1:size(Data,2)
master_data = [master_data, coupling_func(Data(master_dim,i))];
end
slave_data = Data(slave_dim,:);

[Priors, Mu, Sigma] = EM_init_kmeans([master_data;slave_data],K); %EM


%GMM parameters
structGMM=[];
structGMM.Priors = Priors;
structGMM.Mu = Mu;
structGMM.Sigma=Sigma;

% Additional information required at runtime
structGMM.master_dim = master_dim;
structGMM.slave_dim = slave_dim;
structGMM.cplfunc = coupling_func;

% tar = GMR(structGMM.priors, structGMM.mu, structGMM.sigma, 0,1,2:dim);










