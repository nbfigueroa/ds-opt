function [xd b_contour M] = obs_modulation_ellipsoid(x,xd,obs,b_contour,xd_obs)
%
% Obstacle avoidance module: Version 1.1, issued on March 26, 2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Copyright (c) 2011 S. Mohammad Khansari-Zadeh, LASA Lab, EPFL,    %%%
%%%          CH-1015 Lausanne, Switzerland, http://lasa.epfl.ch         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the sufficient modulation
% due to the presence of obstacle(s) so that the generated trajectories
% do not penetrate into the obstacles. 
%
% To avoid struggling with the MATLAB symbolic toolbox, and for the sake of
% simplicity, the provided source code can only be used to avoid obstacles
% in the form of:
%           \Gamma(xt):   \sum_{i=1}^d (xt_i/a_i)^(2p_i) = 1
% For other forms of obstacle shape, it is necessary to modify this file to
% properly compute the matrix of eigenvalues and eigenvectors of the
% dynamic modulation matrix. 
%
% The function is called using:
%       [xd b_contour M] = obs_modulation_ellipsoid(x,xd,obs,b_contour,xd_obs)
%
%
% Inputs -----------------------------------------------------------------
%
%   o x:         d x 1 column vector corresponding to the current robot state
%                probabilities of the K GMM components.
%
%   o xd:        d x 1 column vector corresponding to the derivative of the robot state
%
%   o obs:       A cell array of length N, containing the definition of the N
%                presented obstacles. For example, obs{i} provides the
%                properties of the i-th obstacle. Each obs{i} is a
%                structure variable, and should contains at least the
%                following properties: 
%           - .a:    The ellipsoid length scale (e.g. [1;1] for a 2D circle.
%           - .p:    The ellipsoid power terms (e.g. [1;1] for a 2D circle.
%           - .x0:   The ellipsoid reference point (e.g. [0;0]).
%           - .sf:   The safety factor (optional, default = 1).
%           - .tailEffect (optional, default = true): If it is set to true,
%                    the obstacle modifies the motion even when the it is
%                    moving away from the obstacle. 
%           - .rho (optional, default = 1.0): Sets the reactivity
%                    parameter. The larger the rho, the earlier the robot
%                    responds to the presence of the obstacle. 
%           - .th_r: The obstacle rotation with respect to the global
%                    frame of reference (optional, default = 0.0 rad). 
%           - .partition: Defines the partition of the ellipsoid (optional,
%                    default = [-pi pi]) 
%           Please run 'Tutorial_Obstacle_Avoidance.m' for further information
%           on how to use this obstacle avoidance module.
%
%   o b_contour: A boolean indicating whether the algorithm is in the
%                contouring stage or not.
%
%   o xd_obs:    d x 1 column vector defining the obstacle velocity
%
% Outputs ----------------------------------------------------------------
%
%   o xd:        d x 1 column vector corresponding to the modulated
%                robot velocity.
%
%   o b_contour: A boolean indicating whether the algorithm is in the
%                contouring stage or not.
%
%   o M:         d x d matrix representing the dynamic modulation matrix.
% 
% 
% This code is writen based on the following paper:
% 
%     S.M. Khansari Zadeh and A. Billard, "A Dynamical System Approach to
%     Realtime Obstacle Avoidance", Autonomous Robots, 2012
%
%%
N = length(obs); %number of obstacles
d = size(x,1);
Gamma = zeros(1,N);

xd = xd-xd_obs ; %computing the relative velocity with respect to the obstacle

if d==3
    E = zeros(d,d+1,N);
else
    E = zeros(d,d,N);
end
R = zeros(d,d,N);
M = eye(d);

for n=1:N
    % rotating the query point into the obstacle frame of reference
    if isfield(obs{n},'th_r')
        R(:,:,n) = compute_R(d,obs{n}.th_r);
    else
        R(:,:,n) = eye(d);
    end
    x_t = R(:,:,n)'*(x-obs{n}.x0);
    [E(:,:,n) Gamma(n)] = compute_basis_matrix(d,x_t,obs{n});
end

% [tmp, obs_order] = sort(Gamma,'descend');
w = compute_weights(Gamma,N);
obs_order = 1:N;

for n = obs_order;
    if isfield(obs{n},'rho')
        rho = obs{n}.rho;
    else
        rho = 1;
    end
    D = w(n)*([-1;ones(size(E,2)-1,1)]/abs(Gamma(n))^(1/rho));
    
    if isfield(obs{n},'tailEffect') && ~obs{n}.tailEffect && xd'*R(:,:,n)*E(:,1,n)>=0 %the obstacle is already passed, no need to do anything
        D(1) = 0.0;
    end
    
    if D(1) < -1.0
        D(2:end) = 1.0;
        if xd'*R(:,:,n)*E(:,1,n) < 0
            D(1) = -1.0;
        else
            D(1) = -1.0;
        end
    end
    
    M = (R(:,:,n)*E(:,:,n)*diag(D+1)/E(:,:,n)*R(:,:,n)')*M;
end

E(:,:,n) = R(:,:,n)*E(:,:,n); %transforming the basis vector into the global coordinate system

if b_contour==0 && (D(1) < -0.98) && (E(:,1,n)'*xd < 0) && (norm(M*xd)<0.02)
    b_contour = true;
    disp('Contouring started ... ')
end

if b_contour==1
    contour_dir = sum(E(:,obs{n}.extra.ind,n),2); %extra.ind defines the desired eigenvalues to move along it
    contour_dir = contour_dir/norm(contour_dir);
   
    if (contour_dir'*M*xd >0 && norm(M*xd) > 0.05) || (xd'*E(:,1,n)>0)
        b_contour = false;
        xd = M*xd; %velocity modulation
        disp('Contouring stopped.')
    else
        xd = obs{n}.extra.C_Amp*contour_dir; %extra.C_Amp is the desired amplitude of movement along the controuring direction
    end
else
    xd = M*xd; %velocity modulation
end

xd = xd + xd_obs ; %transforming back the velocity into the global coordinate system

function [E Gamma] = compute_basis_matrix(d,x_t,obs)
% For an arbitrary shap, the next two lines are used to find the shape segment
th = atan2(x_t(2),x_t(1));
if isfield(obs,'partition')
    ind = find(th>=(obs.partition(:,1)) & th<=(obs.partition(:,2)),1);
else
    ind = 1;
end
if isfield(obs,'sf')
    a = obs.a(:,ind).*obs.sf;
else
    a = obs.a(:,ind);
end
p = obs.p(:,ind);
Gamma = sum((x_t./a).^(2*p));

nv = (2*p.*(x_t./a).^(2*p - 1))./a; %normal vector of the tangential hyper-plane

%generating E, for a 2D model it simply is: E = [dx [-dx(2);dx(1)]];
E = zeros(d,d);
E(:,1) = nv;
E(1,2:d) = nv(2:d)';
E(2:d,2:d) = -eye(d-1)*nv(1);

if d == 3
    E(:,end+1) = [0;-nv(3);nv(2)];
end


function w = compute_weights(Gamma,N)
w = zeros(1,N);
Gamma(Gamma<1) = 1;
Gamma = Gamma-1;
for i=1:N
    ind = 1:N;
    ind(i) = [];
    w(i) = prod(Gamma(ind)./(Gamma(i)+Gamma(ind)));
end


function R = compute_R(d,th_r)
% rotating the query point into the obstacle frame of reference

if d == 2 
    R = [cos(th_r(1)) -sin(th_r(1));sin(th_r(1)) cos(th_r(1))];
elseif d == 3
    R_x = [ 1, 0, 0; 0, cos(th_r(1)), sin(th_r(1)); 0, -sin(th_r(1)), cos(th_r(1))];
    R_y = [cos(th_r(2)), 0, -sin(th_r(2)); 0, 1, 0; sin(th_r(2)), 0, cos(th_r(2))];
    R_z = [cos(th_r(3)), sin(th_r(3)), 0; -sin(th_r(3)), cos(th_r(3)), 0; 0, 0, 1];
    R = R_x*R_y*R_z;
else %rotation is not yet supported for d > 3
    R = eye(d);
end