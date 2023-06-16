set(0,'format','shortg'); set(0,'formatspacing','compact');
x = 3.;
err = 1;
while (err > 1e-5)
    xold = x;
    x = log(2.5*x)/0.75;
    err = (xold - x)^2 % watch the results
end
x % echo the answer