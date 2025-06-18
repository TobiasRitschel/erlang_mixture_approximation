function [I, T, F] = adaptive_quadrature_trapezoidal_recursion(f, tl, tu, fl, fu, delta, varargin)
% REFERENCE:
% https://www-personal.umich.edu/~mejn/cp/chapters/int.pdf
% Most of the method is described in exercise 5.20.

% Midpoint
tm = 0.5*(tl + tu);

% Evaluate function at midpoint
fm = f(tm, varargin{:});

% Interval width
dt = tu - tl;

% Width of subdivided interval
dt_fine = 0.5*dt;

% Integral approximation
I_coarse = 0.5*(fl + fu)*dt;

% Integrals of subdivided domains
Il_fine = 0.5*(fl + fm)*dt_fine;
Iu_fine = 0.5*(fm + fu)*dt_fine;

% Fine integral approximation
I_fine = Il_fine + Iu_fine;

% Approximate error
%  - see (5.28)
e = abs(I_fine - I_coarse)/3;

% Is the error sufficiently low?
ErrorIsSufficientlyLow = (e < delta*dt_fine);
% ErrorIsSufficientlyLow = (e < delta*dt);

if(ErrorIsSufficientlyLow)
    % Reuse computations to approximate integral using Simpson's rule
    %  - see (5.9)
    I = (fl + 4*fm + fu)*dt/6;

    % Return times
    T = [tl, tu];

    % Return values of antiderivative
    F = [0, I];
else
    % Repeat computations for subintervals
    [Il, Tl, Fl] = adaptive_quadrature_trapezoidal_recursion(f, tl, tm, fl, fm, delta, varargin{:});
    [Iu, Tu, Fu] = adaptive_quadrature_trapezoidal_recursion(f, tm, tu, fm, fu, delta, varargin{:});

    % Add integrals
    I = Il + Iu;

    % Combine times
    T = [Tl, Tu(2:end)];

    % Combine antiderivative function evaluations
    F = [Fl, Fl(end) + Fu(2:end)];
end
end