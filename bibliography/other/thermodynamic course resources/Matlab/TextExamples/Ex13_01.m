function Ex13_01()
% This function is an example that calls UNIFAC

% The unifac calculation requires the following
% files in the current directory (also showing the calling method:
% unifacAij.mat (parameter file)
% [paramsPure{:}] = unifacSetUp(compArray)
% [gamma]=unifac(x, T, paramsPure{:}) 

addpath(genpath('../gammaModels'));
addpath(genpath('../Psat'));

% specify the groups (not necessarily the smallest
% groups for a molecule) When functional groups
% are near ends, they are usually kept with 
% the ends. See text for examples. compArray is the 
% listing of names and 'nu' values. Each column of 
% 'nu' values correspond to a compound.

% use "load('../gammaModels/UnifacVLE/unifacAij.mat')" to see groups available

% get pure component vapor pressures;
id = [6 44]; % isopropanol and water
[names A B C] = AntoineGet(id)
% propionic acid is not in database
names{3} = 'propionic acid';
A(3) = 0;
B(3) = 0;
C(3) = 0;

% isopropanol + water + propionic acid
compArray = { 
    'CH3'       2       0       1;
    'CH2'       0       0       1;
    'CH'        1       0       0;
    'OH'        1       0       0;
    'H2O'       0       1       0;
    'COOH'      0       0       1;
    };


% This routine needs to be called only once for each
% list of components. The pure parameters are passed
% back merged into a single cell array.
[paramsPure] = unifacSetUp(compArray);

% set concentration and temperature.
% Don't use exactly zero. Use a very small number (1e-50) for zero.
x = [0.6854 0.3146 1e-50];
% normalize to assure sum to one.
x = x/sum(x)
TCelsius = 80.37

% The unifac function can be called whenever x or T changes.
% The pure component cell array is passed along with x and T.
[gamma] = unifac(x,TCelsius, paramsPure{:})

Psat= 10.^(A - B./(TCelsius + C));
Pbub = x.*gamma*Psat'
y = x.*gamma.*Psat/Pbub

[TCelsius fval exitflag] = fzero(@bubP,TCelsius)

function [obj]=bubP(TCelsius)
    [gamma] = unifac(x,TCelsius, paramsPure{:});
    Psat= 10.^(A - B./(TCelsius + C));
    Pbub = x.*gamma*Psat';
    obj = Pbub - 760;
end
end

% ver 1.01 6/5/12 improve documentation