function th = identify_domain_bisection_adaptive_quadrature_recursion(alpha, epsilon, tl, tu, Il, tol, tol_quadrature, varargin)
% Here, we assume that tl and tu are obtained as the finest grid provided
% by the adaptive quadrature rule. Consequently, all subsequent calls to
% the adaptive quadrature routine should not further subdivide the domain.
% However, note also that, as we subdivide the domain as part of the
% bisection, the approximation becomes more accurate than in other parts of
% the approximation.

% Compute midpoint
tm = 0.5*(tu + tl);

% Evaluate cumulative distribution function at the midpoint
Im = adaptive_quadrature_trapezoidal(alpha, tl, tu, tol_quadrature, varargin{:});

% Approximation of integral at midpoint
I = Il + Im;

% Check for convergence
Converged = abs(1 - I - epsilon) < tol;

% Return if the algorithm has converged
if(Converged), th = tm; return; end

if(1 - I < epsilon)
    % Repeat for the lower interval
    th = identify_domain_bisection_adaptive_quadrature_recursion(alpha, epsilon, tl, tm, Il, tol, tol_quadrature, varargin{:});
else
    % Repeat for the upper interval
    th = identify_domain_bisection_adaptive_quadrature_recursion(alpha, epsilon, tm, tu, Im, tol, tol_quadrature, varargin{:});
end
end