function alpha = mixed_folded_normal_kernel_velocity_profile(t, mu0, sigma0, Ns)
% Weights
w = 1/Ns*ones(Ns, 1);

% Scale parameters
sigma = 1.5.^(0:Ns-1)'*sigma0;
mu    = mu0 + [0; cumsum(sigma(1:Ns-1))];

% Subkernels
alpha_i = 1./sqrt(2*pi*sigma.^2) ...
    .*(exp(-(t - mu).^2./(2*sigma.^2)) + exp(-(t + mu).^2./(2*sigma.^2)));

% Kernel
alpha = sum(w.*alpha_i, 1);