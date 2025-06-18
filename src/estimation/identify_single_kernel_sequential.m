function [c, a, M] = identify_single_kernel_sequential(t, M0, dM, alpha, opts, varargin)
% Evaluate true kernel
alphameas = alpha(t, varargin{:});

% Overwrite specific fmincon settings
fmincon_opts = optimoptions(opts.fmincon_opts,                               ...
    'SpecifyObjectiveGradient', true,                                        ...
    'HessianFcn',               @least_squares_kernel_objective_hessian_mex, ...
    'CheckGradients',           false);

% Stopping tolerance
objective_function_tolerance = 1e-10;
if(~isempty(opts) && isfield(opts, 'objective_function_tolerance'))
    objective_function_tolerance = opts.objective_function_tolerance;
end

%% Estimate kernel parameters
% Inequality constraint system matrix and right-hand side
A = [];
B = [];

% Nonlinear constraint function
nonlcon = [];

% Initialize
M = M0;

% Initial guess of coefficients
c0 = ones(1, M+1)/(M+1);
if(~isempty(opts) && isfield(opts, 'c0')), c0(:) = opts.c0; end

% Initial guess rate parameter
a0 = 20;
if(~isempty(opts) && isfield(opts, 'a0')), a0 = opts.a0; end

% Collect initial guesses
theta0 = [c0, a0];

Converged = false;
while(~Converged)
    % Equality constraint system matrix and right-hand side
    Aeq = [ones(1, M+1), 0];
    Beq = 1;

    % Upper and lower bounds
    ub = [ ones(M+1, 1); inf];
    lb = [zeros(M+1, 1);   1];

    % Estimate the interpolation coefficients
    [thetaest, obj, exitflag, output, lambda, grad, hessian] = ...
        fmincon(@least_squares_kernel_objective_mex, ...
        theta0, A, B, Aeq, Beq, lb, ub, nonlcon, fmincon_opts, ...
        t, alphameas, opts); %#ok

    % Parameters in approximate (mixed Erlang) kernel
    c = thetaest(1:end-1);
    a = thetaest(end);

    % Check for convergence
    Converged = (obj < objective_function_tolerance) && (exitflag == 1);

    if(Converged)
        % Break the loop
        break;
    else
        % Increment order
        M = M + dM;

        % Copy initial guess
        theta0 = [c, 1e-8*ones(1, dM), a];
    end
end