% Lactic Acid oligomers, new example 17.10

% part (a)
[x fval exit] = fzero(@(x)(0.2023*(0.555 - 2*x)^2 - x*(2.775 + x)),0)

% part b
% using fsolve from optimization toolbox
[y fval exit] = fsolve(@(x)[0.2023*(0.888 - 2*x(1) - x(2))^2 - (x(1)-x(2))*(1.110 + x(1) + x(2));
0.2023*(x(1)-x(2))*(0.888 - 2*x(1) - x(2)) - x(2)*(1.110 + x(1) + x(2))],[0 0])

% using fminsearch in standard package
[y fval exit] = fminsearch(@(x)(0.2023*(0.888 - 2*x(1) - x(2))^2 ...
    - (x(1)-x(2))*(1.110 + x(1) + x(2)))^2 + ...
    (0.2023*(x(1)-x(2))*(0.888 - 2*x(1) - x(2)) ...
    - x(2)*(1.110 + x(1) + x(2)))^2,[0 0])