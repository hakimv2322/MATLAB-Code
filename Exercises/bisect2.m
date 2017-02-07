function [x, out] = bisect2(FunFcnIn, Intv, params)
% 
%  To find a solution to f(x) = 0 given the continuous function
%  f on the interval [a,b], where f(a) and f(b) have
%  opposite signs:
% 
%  INPUT:  function f(x) defined by function handle FunFcnIn, 
%          interval [a,b]= [Intv.a, Intv.b]
%          tolerence params.tol, max # of iterations = params.MaxIt
%
% PARAMETERS INPUT:
%  FunFcnIn(function handle): function to find the root of
%  Intv(struct): [Intv.a, Intv.b] is the interval containing the root.
%  params(struct): params.tol is the tolerance, params.MaxIt is the
%  maximum iterations to perform.
%
%  OUTPUT: root x, and data structure out. 
%          The success flag out.flg, is 0 for successful
%          execution and non-zero otherwise. out.it is the number
%          of iterations to reach within tolerance.
%
% Written by Ming Gu for Math128A, Fall 2010

TOL = params.tol;
NO  = params.MaxIt;
[FunFcn,msg] = fcnchk(FunFcnIn,0);
if ~isempty(msg)
    error('InvalidFUN',msg);
end
a    = Intv.a;
b    = Intv.b;
if (a > b)
   a = Intv.b;
   b = Intv.a;
end
fa   = sign(FunFcn(a));
fb   = sign(FunFcn(b));
if (fa*fb >0)
    error('Initial Interval may not contain root',msg);
end
if a==b
    error('Initial values for a and b must not equal',msg);
end

It = 0;
out.x =[a;b];
out.f =[FunFcn(a);FunFcn(b)];
while (It <= NO)
   c = (a+b)/2;
   out.x = [out.x;c];
   out.f =[out.f;FunFcn(c)];
   fc = sign(FunFcn(c));
   if (fc ==0)
      x = c;
      out.flg = 0;
      out.it  = It;
      return;
   end
   if (fc * fa < 0)
      b = c;
   else
      a = c;
   end
   if (abs(b-a)<=TOL)
      x = (a+b)/2;
      out.flg = 0;
      out.it  = It;
      return;
   end
   It = It + 1;
end
out.flg =1;
out.it = NO;
x = (a+b)/2;
