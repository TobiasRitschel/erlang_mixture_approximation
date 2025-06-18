function [I, T, F] = adaptive_quadrature_trapezoidal(f, t0, tf, tol, varargin)
% REFERENCE:
% https://www-personal.umich.edu/~mejn/cp/chapters/int.pdf
% Most of the method is described in exercise 5.20.

% Width of interval
dt = tf - t0;

% Is the interval trivial?
IntervalIsTrivial = (dt < tol);

if(IntervalIsTrivial)
    I = 0;
else
    % Tolerance per time unit
    delta = tol/dt;
    
    % Evaluate function
    f0 = f(t0);
    ff = f(tf);
    
    % Call recursion
    [I, T, F] = adaptive_quadrature_trapezoidal_recursion(f, t0, tf, f0, ff, delta, varargin{:});
end
end