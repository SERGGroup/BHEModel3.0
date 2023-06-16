function fzeroExample
% this is a 1-D solver example
set(0,'format','shortg'); set(0,'formatspacing','compact');
xguess = [0 1]; %initial interval, try [0 1] and [-3 -1]
[x,fval,exitflag] = fzero(@func, xguess) %call fzero
%successful results if exitflag = 1

function obj = func(x)
   obj = x^2 + 2*x - 1; %definition of objective
end
end