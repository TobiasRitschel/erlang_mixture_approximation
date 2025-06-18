function th = identify_domain_bisection_adaptive_quadrature(alpha, epsilon, opts, varargin)
% Function tolerance
tol = 1e-14;
if(~isempty(opts) && isfield(opts, 'bisection_tol')), tol = opts.bisection_tol; end

% Tolerance for quadrature
tol_quadrature = 1e-1*tol;

% Left end point
tl = 0;

% Right end point
tu = 1;

% Interval length
dt = tu - tl;

% Initialize integral
Il = 0;

while(true)
    % Approximate cumulative distribution function
    [Iu, Tu, Fu] = adaptive_quadrature_trapezoidal(alpha, tl, tu, tol_quadrature, varargin{:});

    % Total integral
    I = Il + Iu;
    
    % Check for convergence
    Converged = 1 - I < epsilon;

    % Update domain and integral
    if(Converged), break; else, tl = tu; tu = tu + dt; Il = I; end
end

% Identify points bracketing the domain boundary
idxu = find(1 - Il - epsilon < Fu, 1, 'First');
idxl = idxu - 1;

% Identify times bracketing the domain boundary
tl = Tu(idxl);
tu = Tu(idxu);

% Integral value at the lower bracketing point in time
Il = Il + Fu(idxl);

% Solve for the identified domain using bisection
th = identify_domain_bisection_adaptive_quadrature_recursion(alpha, epsilon, tl, tu, Il, tol, tol_quadrature, varargin{:});
end