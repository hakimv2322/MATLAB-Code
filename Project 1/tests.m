function relErrors = tests(C)
% This is a test for the homemade quartic function.
% The quartic function finds all complex roots of
% x^4 + a*x^3 + b*x^2 + c*x - 1 = 0.
% Input is C = [a, b, c].
% It outputs the relative errors of the quartic
% function with the MATLAB roots function.
% (Rel. Error) = |exp - theor|/|theor|

format long
theor = roots([1, C(1), C(2), C(3), -1]);
%Alter the line above to test the local cubic function.
exp = quartic24269728(C);
leng = length(theor);

for (it1 = 1:leng)
    temp = zeros(1, leng);
    for (it2 = 1:leng)
        temp(it2) = abs(exp(it1) - theor(it2))/abs(theor(it2));
    end
    error = min(temp);
    fprintf('Root %d relative error:\n', it1);
    disp(error);
end

end