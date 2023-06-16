function [gamma1, gamma2] = Marg1P(A12,x1, x2)
% calculates gammas for binary using 1 Param Margules
% A12 is a scalar parameter.
% x1 and x2 are mole fractions.
% This function can be used with scalar or vector x1, x2.
% When multiple compositions are passed, the elements of
% x1(i) and x2(i) should be pairs. 
% The gammas have the same dimensions as the x1 and x2 and
% the gamma1(i) and gamma2(i) are the pair that matches 
% x1(i) and x2(i).
gamma1 = exp(A12*(x2).^2);
gamma2 = exp(A12*(x1).^2);
return

% ver. 1.01 4/19/13 improve documentation.