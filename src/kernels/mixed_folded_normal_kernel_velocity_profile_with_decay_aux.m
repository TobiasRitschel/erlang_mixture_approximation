function [w, mu, sigma, normalization_constant] = ...
    mixed_folded_normal_kernel_velocity_profile_with_decay_aux ...
    (mu0, sigma0, lambda, Ns)
% Weights
w = 1/Ns*ones(Ns, 1);

% Scale and location parameters
sigma = 1.5.^(0:Ns-1)'*sigma0;
mu    = mu0 + [0; cumsum(sigma(1:Ns-1))];

% Normalization constant
normalization_constant = 0.5*exp(0.5*sigma.^2*lambda^2) ...
    .*(exp( lambda*mu).*erfc((lambda*sigma.^2 + mu)./(sqrt(2)*sigma)) ...
    +  exp(-lambda*mu).*erfc((lambda*sigma.^2 - mu)./(sqrt(2)*sigma)));