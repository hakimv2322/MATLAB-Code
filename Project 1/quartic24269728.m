function rts = quartic24269728(C)

% This function solves for all complex roots of
% the quartic:
% x^4 + a*x^3 + b*x^2 + c*x - 1 = 0
% (notice this polynomial always has a root)
%
% The input vector is C = [a, b, c], where
% a, b, c are all real coefficients.
%
% Roots of multiplicity m are in the output m times.
%
% The strategy is to use Newton's Method starting
% at the "edge" of the polynomial, where there
% are no roots either to the left or to the
% right of the starting point.

format long
quad([1.527652, 2.221, 0.0000001])

end

function rts = quad(C)
% This solves for the complex roots of
% the quadratic
% a*x^2 + b*x + c = 0
%
% The input vector is C = [a, b, c], where
% a, b, c are all real coefficients. a is nonzero.
% Double roots are output twice.

a = C(1); b = C(2); c = C(3);
delta = sqrt(b^2 - 4*a*c);
if (b >= 0)
    x1 = (-b - delta)/(2*a);
    x2 = c/(a*x1);
    rts = [x1, x2];
else
    x1 = (-b + delta)/(2*a);
    x2 = c/(a*x1);
    rts = [x1, x2];

end
end

