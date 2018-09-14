function [x_dot] = lin_ds(b,x,type,alpha,varargin)
if nargin > 4
    xi_0 = varargin{1};
else
    xi_0 = [-2;0.2];
end

switch type
    case 1
    % Linear DS converging to a target
    A = -[10 0;0 10];
    
    case 2
    % Linear DS diverging from a target
    A = [10 0;0 10];
    
    case 3
    % Linear DS on x
    A = [-10 0;0 0];
    x_dot = A*(x - repmat([1 0]', [1 size(x,2)]));
    
    case 4
    % Linear DS on y
    A = [0 0;0 -10];
    x_dot = A*(x + repmat([0 1]', [1 size(x,2)]));
    
    case 5
    % Linear DS diagonal
    A = [10 0;0 -10];
    
    case 6
    % Linear DS diagonal target 1
    A = [0 -1;1 -2];
    
    case 7
    % Linear DS diagonal target 2
    A = [0 1;-1 -2];    
    
    case 8
    % Linear DS diagonal target 2
    A = [0 -1; 1 -2];    
    
    case 9
    % Linear DS diagonal target 2
    A = [0 1; -1 -2];    
    
    case 10
    % Linear DS sim. ref lin. traj. with target as eig -> lambda -1
     y1 = 1;
     y2 = -xi_0(1)/xi_0(2);
     y = [y1;y2];    
     Q = [y./norm(y),xi_0./norm(xi_0)];
     L = [-10 0 ; 0 -1];
     A = Q*L*Q';               
     
    case 11
     % Linear DS converging to a target
     A1 = -[10 0;0 10];
     
     % Linear DS sim. ref lin. traj. with target as eig -> lambda -1
     y1 = 1;
     y2 = -xi_0(1)/xi_0(2);
     y = [y1;y2];    
     Q = [y./norm(y),xi_0./norm(xi_0)];
     L = [-10 0 ; 0 -1];
     A2 = Q*L*Q';               
     A = alpha*A1 + (1-alpha)*A2;
end

x_dot = A*x;

% if (type>4) || (type < 3) 
%     if isempty(b)
%         x_dot = A*x;
%     else
%         x_dot = A*x + repmat(b,[1 length(x)]);
%     end
% end

end