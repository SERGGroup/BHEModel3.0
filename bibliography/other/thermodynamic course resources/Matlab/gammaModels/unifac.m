function [gamma] = unifac(x, TCelsius, r, q, R, Q, nComp, nGroups, compNArray, X, THETA, aij)
% This function calcualtes UNIFAC activity coefficients for binary and multicomponent mixtures.
% Only one composition can be evaluated with each function call.
% x - mole fraction row vector, with n elements for an n-compoent mixture.
% TCelsius - system temperature in Celsius.
% gamma - row vector of activity coefficients in the same order as x.
% The common way this function is called is
%  [gamma] = unifac(x,TCelsius, paramsPure{:});
% When distributed, the pure parameters
% r, q, R, Q, nComp, nGroups, compNArray, X, THETA, aij
% are loaded into a single cell array 'paramsPure' by
% using a call to unifacSetup.m before this funciton is called.
% Only three arguments are passed when the pure parameters
% are combined into the cell array.
% The cell array is automatically 'dealt' to the 
% corresponding variable list when the function runs.

% variables are the same as textbook with the following additions:
% compNArray - the numerical portion of the compArray
%              listing the number of groups in each compound.
%              This is a table of 'nu' values for each group.

T = TCelsius + 273.15;
% make sure no x value is identically zero
x = x + 1e-50;

% check that x is the correct size for the pure component list
sizex = size(x);
if sizex(1) ~= 1
    fprintf('Composition vector x should be a single row vector, but has %d rows. \nErrors will result.\n',sizex(1))
end
if sizex(2) ~= size(r,2)
    fprintf('Composition vector x has %d columns, but compArray indicates %d components. \nErros will result.\n',sizex(2),size(r,2))
end

% now calculate the sumoi(nu.j*x.i)
% The 'nu' values are in compNArray
sumoi = compNArray*x';

% determine the mol fractions of groups in the solution
% The result is a column vector with mole fractions
% of each group in the same order listed in compNArray
Xmix = sumoi/sum(sumoi)';

% Append Xmix to the end of X table from pure components
% See unifacSetUp for description of X.
X = [X Xmix];

% PSI matrix
PSI = exp(-aij/T);

% calculate group surface area fractions of mixture
% THETAmix is a vector with the surface fractions for each group.
QXmix = Q'.*Xmix;
THETAmix = QXmix/sum(QXmix)';

% append THETAmix to THETA of the pure components.
% See unifacSetUp.m for description of THETA
THETA = [THETA  THETAmix];

% Activity formulas require sums of THETA*PSI
sumoiTHETAiPSIij = PSI'*THETA;
THETAjosumoiTHETAiPSIij = THETA./sumoiTHETAiPSIij;
sumRatio = PSI*THETAjosumoiTHETAiPSIij;

% determine the lnGAMMA for the pure components and mixture
% mixture will be in the last column.
% The kron function is combined with a ones matrix to build
% copies of the Q vector before matix dot multiplication.
LnGAMMA = kron(Q',ones(1,nComp+1)).*(1-log(sumoiTHETAiPSIij) - sumRatio);

%set up storage for residual gammas
lngammaRes= zeros(1,nComp);
% calculate residual gammas, using the last row (mixture) relative
% to pure components. The 'nu' values are in compNArray.
for i = 1:nComp 
    lngammaRes(i) = compNArray(:,i)'*(LnGAMMA(:,nComp+1)-LnGAMMA(:,i));
end

% now calculate the configurational part
theta = x.*q/(x*q');
phi = x.*r/(x*r');
lngammaC = log(phi./x) + (1 - phi./x) - 5 * q .* (log(phi./theta)+ ( 1 - phi./theta));

% the vector of activity coefficients by combining combinitorial and
% residual parts.
gamma = exp(lngammaC + lngammaRes);
end

% ver 1.02 4/19/13 added documentation.
% ver 1.01 added dimension checking for input arrays.