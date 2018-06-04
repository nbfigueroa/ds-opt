function [A_c, b_c, P] = optimize_lpv_ds_from_data(Data, attractor, ctr_type, gmm, varargin)

% Positions and Velocity Trajectories
Xi_ref = Data(1:2,:);
Xi_ref_dot = Data(3:4,:);
[N,M] = size(Xi_ref_dot);

% Define Optimization Variables
sdp_options = []; Constraints = [];
epsilon = 0.001;
init_cvx = 0;

% Define DS Variables
K = length(gmm.Priors);
A_c = zeros(N,N,K);
b_c = zeros(N,K);

% Define solver for type of constraints
switch ctr_type
    case 0
        % 'sedumi': semidefinite programming solver for convex problems
        sdp_options = sdpsettings('solver','sedumi','verbose', 1);
    
    case 1
        % 'penlab': Nonlinear semidefinite programming solver
        sdp_options = sdpsettings('solver','penlab','verbose', 1,'usex0',1);
        P_var = sdpvar(N, N, 'symmetric','real');
        Constraints = [Constraints, P_var >  eye(N,N)];
        assign(P_var,eye(N));
        init_cvx = varargin{2};
        
    case 2
        % 'penlab': Nonlinear semidefinite programming solver
        sdp_options = sdpsettings('solver','penlab','verbose', 1,'usex0',1);
        P = varargin{1};
        init_cvx = varargin{2};
end

if init_cvx
    % Solve Problem with Convex constraints first to get A's
    fprintf('Solving Optimization Problem with Convex Constraints for Non-Convex Initialization...\n');
    [A0, b0] = optimize_lpv_ds_from_data(Data, attractor, 0, gmm);
end

% Posterior Probabilities per local model
h_k = posterior_probs_gmm(Xi_ref,gmm,'norm');


% Define Constraints and Assign Initial Values
for k = 1:K    
    A_vars{k} = sdpvar(N, N, 'full','real');       
    b_vars{k} = sdpvar(N, 1, 'full');
    Q_vars{k} = sdpvar(N, N,'symmetric','real');       
       
    % Assign Initial Parameters
    if init_cvx       
        assign(A_vars{k},A0(:,:,k));
        assign(b_vars{k},b0(:,k));        
    else
        assign(A_vars{k},-eye(N));
        assign(b_vars{k},-eye(N)*attractor);
    end
    
    % Define Constraints
    switch ctr_type
        case 0 %: convex
            Constraints = [Constraints transpose(A_vars{k}) + A_vars{k} <= -epsilon*eye(N,N)]; 
            Constraints = [Constraints b_vars{k} == -A_vars{k}*attractor];
        
        case 1 %: non-convex, unknown P                                                        
            Constraints = [Constraints, transpose(A_vars{k})*P_var + P_var*A_vars{k} <= -epsilon*eye(N)];
            Constraints = [Constraints, b_vars{k} == -A_vars{k}*attractor];
            

         
        case 2 %: non-convex with given P
                        
            % Option 1: Stricter Constraint                      
%             Constraints = [Constraints, transpose(A_vars{k})*P + P*A_vars{k} <= -epsilon*eye(N)];                        
            
            % Option 2: Less Strict and converges faster most of the times                      
            Constraints = [Constraints, transpose(A_vars{k})*P + P*A_vars{k} == Q_vars{k}];                        
            Constraints = [Constraints, Q_vars{k} <= -epsilon*eye(N)];                        
            Constraints = [Constraints, b_vars{k} == -A_vars{k}*attractor];                                 
            assign(Q_vars{k},-eye(N));
            
    end
end

% Calculate our estimated velocities caused by each local behavior
Xi_d_dot_c_raw = sdpvar(N,M,K, 'full');%zeros(size(Qd));
for k = 1:K
    h_K = repmat(h_k(k,:),[N 1]);
    f_k = A_vars{k}*Xi_ref + repmat(b_vars{k},[1 M]);
    Xi_d_dot_c_raw(:,:,k) = h_K.*f_k;
end

% Sum each of the local behaviors to generate the overall behavior at
% each point
Xi_d_dot = sdpvar(N, M, 'full');
Xi_d_dot = reshape(sum(Xi_d_dot_c_raw,3),[N M]);

% Then calculate the difference between approximated velocities
% and the demonstated ones for A
Xi_dot_error = Xi_d_dot - Xi_ref_dot;

% Defining Objective Function depending on constraints
if ctr_type == 0
    Xi_dot_total_error = sdpvar(1,1); Xi_dot_total_error(1,1) = 0;
    for m = 1:M
        Xi_dot_total_error = Xi_dot_total_error + norm(Xi_dot_error(:, m));
    end
    Objective = Xi_dot_total_error;
else
    Aux_var     = sdpvar(N,length(Xi_dot_error));
    Objective   = sum((sum(Aux_var.^2)));
    Constraints = [Constraints, Aux_var == Xi_dot_error];
end

% Solve optimization problem
sol = optimize(Constraints, Objective, sdp_options)
if sol.problem ~= 0
    yalmiperror(sol.problem);
end

for k = 1:K
    A_c(:,:,k) = value(A_vars{k});
    b_c(:,k)   = value(b_vars{k});
end

switch ctr_type 
    case 0 
        P = eye(N);
    case 1
        P = value(P_var);
end

sol.info
check(Constraints)
fprintf('Total error: %2.2f\nComputation Time: %2.2f\n', value(Objective),sol.solvertime);


%%%% FOR DEBUGGING: Check Negative-Definite Constraint %%%%
if ctr_type == 0
    constr_violations = zeros(1,K);
    for k=1:K
        A_t = A_c(:,:,k) + A_c(:,:,k)';
        constr_violations(1,k) = sum(eig(A_t) > 0); % sufficient
    end
    % Check Constraint Violation
    if sum(constr_violations) > 0
        warning(sprintf('Strict System Matrix Constraints are NOT met..'))
    else
        fprintf('All Sufficient System Matrix Constraints are met..\n')
    end
else
    suff_constr_violations = zeros(1,K);
    for k=1:K
        Pg_A =  A_c(:,:,k)'*P + P*A_c(:,:,k);
        suff_constr_violations(1,k) = sum(eig(Pg_A + Pg_A') > 0); % strict
    end
    % Check Constraint Violation
    if sum(suff_constr_violations) > 0
        warning(sprintf('Sufficient System Matrix Constraints are NOT met..'))
    else
        fprintf('All Sufficient System Matrix Constraints are met..\n')
    end
end

end