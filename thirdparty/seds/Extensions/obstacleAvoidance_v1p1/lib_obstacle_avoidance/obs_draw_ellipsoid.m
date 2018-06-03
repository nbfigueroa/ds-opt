function [x_obs x_obs_sf] = obs_draw_ellipsoid(obs,ns)
%
% Obstacle avoidance module: Version 1.1, issued on March 26, 2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Copyright (c) 2011 S. Mohammad Khansari-Zadeh, LASA Lab, EPFL,    %%%
%%%          CH-1015 Lausanne, Switzerland, http://lasa.epfl.ch         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the boundaries of the obstacle, and also the
% border of the safety region. This function is only useful for
% illustration of 2D and 3D obstacles, and is not required for the obstacle
% avoidance.
%
% To avoid struggling with the MATLAB symbolic toolbox, and for the sake of
% simplicity, the provided source code can only be used to avoid obstacles
% in the form of:
%           \Gamma(xt):   \sum_{i=1}^d (xt_i/a_i)^(2p_i) = 1
% For other forms of obstacle shape, it is necessary to modify this file to
% properly compute the border of the obstacle.
%
% The function is called using:
%       [x_obs x_obs_sf] = obs_draw_ellipsoid(obs,ns)
%
%
% Inputs -----------------------------------------------------------------
%
%   o obs:  A cell array of length N, containing the definition of the N
%           presented obstacles. For example, obs{i} provides the
%           properties of the i-th obstacle. Each obs{i} is a
%           structure variable, and should contains at least the
%           following properties: 
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
%           Please run 'Tutorial_Obstacle_Format.m' for further information
%           on how to define obstacles.
%
%   o ns:    Number of segments used to compute the boundary of the
%            obstacle. For 2D obstacles, ns is a positive scalar defining
%            the number of datapoint on the perimeter of the obstacle. For
%            3D obstacles, ns is a 2 x 1 vector. The first component is the
%            number of points in x-y plane, and the second component
%            corresponds to the number of point in semi-half of the y-z
%            plane.
%
% Outputs ----------------------------------------------------------------
%
%   o x_obs:     A d x tp x N matrix containing the value of boundary
%                points of the all obstacles. For 2D obstacles tp = ns,
%                while for 3D obstacles tp = ns(1)*ns(2). 
%                You can use the following codes to illustrate the n-th
%                obstacle:
%
%                For 2D obstacles: 
%                     patch(x_obs(1,:,n),x_obs(2,:,n),[0.6 1 0.6])
%
%                For 3D obstacles:
%                     h = surf(reshape(x_obs(1,:,n),ns(2),ns(1)),...
%                              reshape(x_obs(2,:,n),ns(2),ns(1)),...
%                              reshape(x_obs(3,:,n),ns(2),ns(1))); 
%
%   o x_obs_sf: Same as x_obs, with the difference that it corresponds to
%               the boundary of the safety region.
%
% This code is writen as an accompanying file for the following paper:
% 
%     S.M. Khansari Zadeh and A. Billard, "A Dynamical System Approach to
%     Realtime Obstacle Avoidance", Autonomous Robots, 2012
%
%%
d = size(obs{1}.x0,1);
if d > 3 || d < 2
    disp('Error: obs_draw_ellipsoid only supports 2D or 3D obstacle')
    x_obs = [];
    x_obs_sf = [];
    return
end

if d == 2
    theta = linspace(-pi,pi,ns(1));
else
    [theta phi] = meshgrid(linspace(-pi,pi,ns(1)),linspace(-pi/2,pi/2,ns(2))); %
    theta = theta(:)';
    phi = phi(:)';
end

np = prod(ns);

N = length(obs); %number of obstacles
x_obs = zeros(d,np,N);

if nargout > 1;
    x_obs_sf = zeros(d,np,N);
end
    
for n=1:N
    clear ind
    % rotating the query point into the obstacle frame of reference
    if isfield(obs{n},'th_r')
        if size(obs{n}.th_r(:),1) == d && size(obs{n}.th_r(:),2) == d
            R = obs{n}.th_r;
        else
            R = compute_R(d,obs{n}.th_r);
        end
    else
        R = eye(d);
    end

    % For an arbitrary shap, the next two lines are used to find the shape segment
    if isfield(obs{n},'partition')
        for i=1:size(obs{n}.partition,1)
            ind(i,:) = theta>=(obs{n}.partition(i,1)) & theta<=(obs{n}.partition(i,2));
        end
        [i ind]=max(ind);
    else
        ind = 1;
    end
    a = obs{n}.a(:,ind);
    p = obs{n}.p(:,ind);

    if d == 2
        x_obs(1,:,n) = a(1,:).*cos(theta);
        x_obs(2,:,n) = a(2,:).*sign(theta).*(1 - cos(theta).^(2.*p(1,:))).^(1./(2.*p(2,:)));
    else
        x_obs(1,:,n) = a(1,:).*cos(phi).*cos(theta);
        x_obs(2,:,n) = a(2,:).*sign(theta).*cos(phi).*(1 - 0.^(2.*p(3,:)) - cos(theta).^(2.*p(1,:))).^(1./(2.*p(2,:)));
        x_obs(3,:,n) = a(3,:).*sign(phi).*(1 - (sign(theta).*cos(phi).*(1 - 0.^(2.*p(3,:)) - cos(theta).^(2.*p(1,:))).^(1./(2.*p(2,:)))).^(2.*p(2,:)) - (cos(phi).*cos(theta)).^(2.*p(1,:))).^(1./(2.*p(3,:)));
    end

    if nargout > 1;
        if isfield(obs{n},'sf')
            if length(obs{n}.sf) == 1
                x_obs_sf(:,:,n) = R*(x_obs(:,:,n).*obs{n}.sf) + repmat(obs{n}.x0,1,np);
            else
                x_obs_sf(:,:,n) = R*(x_obs(:,:,n).*repmat(obs{n}.sf,1,np)) + repmat(obs{n}.x0,1,np);
            end
        else
            x_obs_sf(:,:,n) = R*x_obs(:,:,n) + repmat(obs{n}.x0,1,np);
        end
    end
    
    x_obs(:,:,n) = R*x_obs(:,:,n) + repmat(obs{n}.x0,1,np);
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