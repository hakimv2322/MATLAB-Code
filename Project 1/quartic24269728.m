function rts = quartic24269728(C)

% This function solves for all complex rts of
% the quartic:
% x^4 + b*x^3 + c*x^2 + d*x - 1 = 0
% (notice this polynomial always has a rt)
%
% The input vector is C = [b, c, d], where
% b, c, d are all real coefficients.
%
% rts of multiplicity m are in the output m times.
%
% The strategy is to use Newton's Method starting
% at the "edge" of the polynomial, where there
% are no rts of the derivative either to the left
% or to the right of the starting point.

format long
%rts = double(round(cubic(C), 18, 'significant').');
% Use the above to test the local function cubic()

b = double(C(1)); c = double(C(2)); d = double(C(3));
derivRts = cubic([4, 3*b, 2*c, d]);
syms x;
f = poly2sym([1, b, c, d, -1], x);
p0 = derivRts(1);

for (it = 2:3)
    if (imag(derivRts(it)) == 0)
        x = derivRts(it);
        val1 = subs(f);
        x = p0;
        val2 = subs(f);
        if (val2 > val1)
            p0 = derivRts(it);
        end
    end
end

temp = sort(derivRts);
if (p0 == temp(1))
    if (p0 <= 0)
        x1 = Newton([1, b, c, d, -1], 1.1*p0 - 0.01);
    else
        x1 = Newton([1, b, c, d, -1], 0.9*p0);
    end
else
    if (p0 >= 0)
        x1 = Newton([1, b, c, d, -1], 1.1*p0 + 0.01);
    else
        x1 = Newton([1, b, c, d, -1], 0.9*p0);
    end
end

Qb = x1 + b;
Qc = x1^2 + x1*b + c;
Qd = 1/x1;

temp = cubic([1, Qb, Qc, Qd]);
x2 = temp(1);
x3 = temp(2);
x4 = temp(3);

rts = [x1, x2, x3, x4].';

end

function rt = Newton(F, p0)
% This function is Newton's Method,
% for polynomials. Output is the single
% resulting real rt.
%
% p0 is the starting point, a real number.
% F is a vector of coefficients; the first
% entry is for the highest degree term.

counter = 0;

syms x;
f = poly2sym(F, x);
Df = diff(f);
p = p0;
relError = 1.0;

while (relError >= 10^-20)
    counter = counter + 1;
% Uncomment the line below to see each iteration.
%    disp(counter);
    if (counter > 10^6)
        disp('Too many iterations--terminated.')
        disp(double(p));
        break
    end
    
    temp = p;
    x = p;
    if (subs(Df) == 0)
        disp('Attempted a divide by zero.')
        disp(double(p));
        break
    end
    p = p - (subs(f))/(subs(Df));
    p = round(double(p), 24, 'significant');
    relError = abs(temp - p)/abs(p);
end

%disp(counter);
% Uncomment the line above to see how many iterations
% were used in this Newton's Method.

rt = p;

end


function rts = quad(C)
% This solves for the complex rts of
% the quadratic
% a*x^2 + b*x + c = 0
%
% The input vector is C = [a, b, c], where
% a, b, c are all real coefficients. a is nonzero.
% double rts are output twice.

a = double(C(1)); b = double(C(2)); c = double(C(3));
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

function rts = cubic(C)
% This solves for all complex rts of the cubic
% a*x^3 + b*x^2 + c*x + d = 0
%
% The input vector is C = [a, b, c, d], where
% a, b, c, d are all real coefficients. a is nonzero.
% double rts are output twice.

a = double(C(1)); b = double(C(2));
c = double(C(3)); d = double(C(4));
derivRts = quad([3*a, 2*b, c]);
if (imag(derivRts(1)) ~= 0)
    x1 = Newton(C, 0.0);
else
    if (derivRts(1) == derivRts(2))
        Da = derivRts(1);
        if (Da*a >= 0)
            if (Da <= 0)
                x1 = Newton(C, 1.1*Da - 0.01);
            else
                x1 = Newton(C, 0.9*Da);
            end
        else
            if (Da >= 0)
                x1 = Newton(C, 1.1*Da + 0.01);
            else
                x1 = Newton(C, 0.9*Da);
            end
        end
    else
        derivRts = sort(derivRts);
        Da = derivRts(1); Db = derivRts(2);
        fOfDa = a*Da^3 + b*Da^2 + c*Da + d;
        fOfDb = a*Db^3 + b*Db^2 + c*Db + d;
        if (fOfDa*fOfDb <= 0)
            if (Da <= 0)
                x1 = Newton(C, 1.1*Da - 0.01);
            else
                x1 = Newton(C, 0.9*Da);
            end
        else
            if (fOfDa*a > 0)
                if (Da <= 0)
                    x1 = Newton(C, 1.1*Da - 0.01);
                else
                    x1 = Newton(C, 0.9*Da);
                end
            else
                if (Db >= 0)
                    x1 = Newton(C, 1.1*Db + 0.01);
                else
                    x1 = Newton(C, 0.9*Db);
                end
            end
        end
    end
end

if (x1 == 0)
    Qa = a; Qb = b; Qc = c;
else
    Qa = a;
    Qb = -(d + x1*c)/(x1^2);
    Qc = -d/x1;
end
temp = quad([Qa, Qb, Qc]);
x2 = temp(1);
x3 = temp(2);

rts = [x1, x2, x3];

end

