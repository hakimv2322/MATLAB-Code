function relErrors = tests(C)
% This is a test for the homemade quartic function.
% The quartic function finds all complex roots of
% x^4 + a*x^3 + b*x^2 + c*x - 1 = 0.
% Input is C = [a, b, c].
% It outputs the relative errors of the quartic
% function with the MATLAB roots function.
% (Rel. Error) = |exp - theor|/|theor|

format long
theor = sort(roots(C));
exp = sort(quartic24269728(C));
leng = length(theor);
if (imag(theor(leng-1)) > 0)
    temp = theor(leng-1);
    theor(leng-1) = theor(leng);
    theor(leng) = temp;
end
if (imag(exp(leng-1)) > 0)
    temp = exp(leng-1);
    exp(leng-1) = exp(leng);
    exp(leng) = temp;
end

for (it = 1:leng)
    fprintf('Root %d relative error:\n', it);
    disp(abs(exp(it) - theor(it))/abs(theor(it)));
end

end