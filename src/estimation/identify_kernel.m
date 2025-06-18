function [c, a, M, th] = identify_kernel(alpha, beta, opts, varargin)

% Number of kernels
nz = numel(alpha);

% Tolerance for main support
epsilon = 1e-14;
if(~isempty(opts) && isfield(opts, 'epsilon' )), epsilon  = opts.epsilon; end

% Number of measurements
N = 100;
if(~isempty(opts) && isfield(opts, 'N' )), N  = opts.N; end

% Initial guess of integer kernel parameter
M0 = 0;
if(~isempty(opts) && isfield(opts, 'M0')), M0 = opts.M0; end

% Increment in integer kernel parameter search
dM = 5;
if(~isempty(opts) && isfield(opts, 'dM')), dM = opts.dM; end

for i = 1:nz
    if(isempty(beta{i}))
        % Identify domain using bisection
        th = identify_domain_bisection_adaptive_quadrature(alpha{i}, epsilon, opts, varargin{:});
    else
        % Identify domain using bisection
        th = identify_domain_bisection(beta{i}, epsilon, opts, varargin{:});
    end

    % Measurement points
    tmeas = linspace(0, th, N+1);

    % Identify kernel for fixed order
    [ci, ai, Mi] = identify_single_kernel_sequential(tmeas, M0, dM, alpha{i}, opts, varargin{:});

    % Store result
    c(i, 1:Mi+1) = ci; %#ok
    a(i, 1     ) = ai; %#ok
    M(i, 1     ) = Mi; %#ok
end