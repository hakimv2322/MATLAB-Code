function rts = quad_solve(coeff)
% Quadratic solver with 5 sig-fig arithmetic
% coeff is an array with three entries.
a = coeff(1);
b = coeff(2);
c = coeff(3);
% ^MATLAB indexing is non-standard.

delta = fl_5(sqrt(fl_5(fl_5(b^2) - fl_5(4*fl_5(a*c)))));
x1 = fl_5((fl_5(-b + delta))/fl_5(2*a));
x2 = fl_5((fl_5(-b - delta))/fl_5(2*a));
rts = [x1, x2];