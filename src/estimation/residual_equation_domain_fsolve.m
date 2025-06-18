function [res, dres] = residual_equation_domain_fsolve(lnt, alpha, beta, epsilon, varargin)
% Transform decision variable
t = exp(lnt);

if(nargout > 1)
    % Evaluate probability density function (i.e., the derivative)
    f = alpha(t, varargin{:});
end

% Evaluate cumulative distribution function
F = beta(t, varargin{:});

% Evaluate residual equation (log-transformed)
res = log(1 - F) - log(epsilon);

if(nargout > 1)
    % Evaluate derivative of time wrt. log-time
    dt_dlnt = t;

    % Evaluate derivative
    dres = -f./(1 - F)*dt_dlnt;
end
end