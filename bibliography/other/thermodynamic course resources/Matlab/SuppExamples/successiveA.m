set(0,'format','shortg'); set(0,'formatspacing','compact');
x = 1;
err = 1;
while (err > 1e-5)
    xold = x;
    x = exp(0.75*x)/2.5;
    err = (xold - x)^2 % watch the results
end
x % echo the answer