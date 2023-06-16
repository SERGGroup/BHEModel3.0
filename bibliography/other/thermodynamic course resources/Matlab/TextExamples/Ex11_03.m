function Ex11_03
% illustration of fitting Pbub, one parameter Margules
addpath(genpath('../../gammaModels')); %add path to gamma models
x1 = 0.6854; x2 = 1-x1;
A = 1;
P1sat = 694; P2sat = 359.9;
[A fval exitflag] = fzero(@calcObj,A);
A
fval
exitflag
  function obj = calcObj(A)
    [gamma1 gamma2] = Marg1P(A,x1,x2);
    Pbub = (x1.*gamma1*P1sat) + (x2.*gamma2*P2sat);
    obj = Pbub-760;
  end
end
