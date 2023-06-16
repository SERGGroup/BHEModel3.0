function [gamma1, gamma2] = Marg2P(A12,A21, x1, x2)
% calculates gammas for binary using 2 Parameter Margules
% x1, x2 can be vector or scalar.
% A12 and A21 are scalar parameters.
% x1 and x2 are mole fractions.
% This function can be used with scalar or vector x1, x2.
% When multiple compositions are passed, the elements of
% x1(i) and x2(i) should be pairs. 
% The gammas have the same dimensions as the x1 and x2 and
% the gamma1(i) and gamma2(i) are the pair that matches 
% x1(i) and x2(i).
gamma1 = exp((x2.^2).*(A12 + 2*(A21 - A12)*x1));
gamma2 = exp((x1.^2).*(A21 + 2*(A12 - A21)*x2));
return

% ver 1.01 4/20/13 improve documentation