function f = frechet_kernel(t, w, mu, sigma)
% Indices of "valid" times
idx = (t > mu);

% Allocate memory
f = zeros(size(t));

% Evaluate shifted and scaled time
s = (t(idx) - mu)/sigma;

% Evaluate kernel for sufficiently large times
f(idx) = w/sigma.*s.^(-1 - w).*exp(-s.^(-w));