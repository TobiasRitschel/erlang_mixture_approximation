function th = identify_domain_bisection(beta, epsilon, opts, varargin)

% Function tolerance
tol = 1e-15;
if(~isempty(opts) && isfield(opts, 'bisection_tol')), tol = opts.bisection_tol; end

% Left end point
tl = 0;

% Right end point
tu = 1;

% Interval length
dt = tu - tl;

while(true)
    % Evaluate cumulative distribution function
    Fu = beta(tu, varargin{:});
    
    % Check for convergence
    Converged = 1 - Fu < epsilon;

    % Update domain
    if(Converged), break; else, tl = tu; tu = tu + dt; end
end

% Solve for the identified domain using bisection
th = identify_domain_bisection_recursion(beta, epsilon, tl, tu, tol, varargin{:});
end