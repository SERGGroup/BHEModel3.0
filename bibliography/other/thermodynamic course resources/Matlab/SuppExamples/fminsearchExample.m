function fminsearchExample
% this is a 2-D solver example
set(0,'format','shortg'); set(0,'formatspacing','compact');
guess = [1 3]; %initial guess, try [-1 3] and [1 3]
[x,fval,exitflag] = fminsearch(@func, guess) %call fzero
% successful results if exitflag = 1
function obj = func(guess)
    x=guess(1);
    y=guess(2);
    obj1 = x^2 + 2*y - 10;
    obj2 = x + y - 4;
    obj = obj1^2 + obj2^2; %definition of objective
end
end