function th = identify_domain_fsolve(beta, alpha, epsilon, opts, varargin)
% Identify the domain of main support of the probability density function
lnth = fsolve(@residual_equation_domain_fsolve, 0, opts, ...
    alpha, beta, epsilon, varargin{:});

% Dominant domain of support
th = exp(lnth);
end