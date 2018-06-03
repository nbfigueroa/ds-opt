function vel = getCDSVelocity(master_dyn, slave_dyn, cds_model, data_point, alpha, beta)
%
% This function calculates the velocity of both the master and slave
% sub-systems, given their dynamic models and the coupling.
%
%
% Inputs -----------------------------------------------------------------
%
%
%   o master_dyn:  SEDS dynamic model of the master sub-system
%
%   o slave_dyn:   SEDS dynamic model of the slave sub-systetm
%
%   o cds_model:   The coupling model between master and slave
%
%   o data_point:  d x 1 array representing the d-dimensional query point 
%                  in state space.
%
%   o alpha:       Scalar parameter for the CDS model
%
%   o beta:        Scalar parameter for the CDS model
% Outputs ----------------------------------------------------------------
%
%   o vel:        d x 1 array representing the calculated velocity at the
%                 query point
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


% Collecting the information about master and slave systems separately
master_var = data_point(cds_model.master_dim);
master_size = length(cds_model.master_dim);
tr_master_var = cds_model.cplfunc(master_var);
slave_size = length(cds_model.slave_dim);
slave_var = data_point(cds_model.slave_dim);

% Velocity of the master sub-system
m_vel = GMR(master_dyn.Priors, master_dyn.Mu, master_dyn.Sigma, master_var,1:master_size, ...
    master_size+1:2*master_size);

% Error of the slave sub-system
s_err = slave_var - GMR(cds_model.Priors, cds_model.Mu, cds_model.Sigma, ...
    tr_master_var, 1:length(tr_master_var), length(tr_master_var)+1:length(tr_master_var)+slave_size);

% Correction to be applied if the coupling function does not pass through
% the origin
st_correc = GMR(cds_model.Priors, cds_model.Mu, cds_model.Sigma, ...
    zeros(size(master_var)), 1:length(tr_master_var), length(tr_master_var)+1:length(tr_master_var)+slave_size);

% Velocity of the slave sub-system
s_vel = alpha*GMR(slave_dyn.Priors, slave_dyn.Mu, slave_dyn.Sigma, beta*(s_err+st_correc), 1:slave_size, slave_size+1:2*slave_size);

%Recollecting both velocities in one vector
vel=zeros(master_size+slave_size,1);
vel(cds_model.master_dim) = m_vel;
vel(cds_model.slave_dim) = s_vel;


