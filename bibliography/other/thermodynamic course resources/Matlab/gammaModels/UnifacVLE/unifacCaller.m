% This script is an example that calls UNIFAC

% The unifac calculation requires the following
% files (also showing the calling method):
% unifacAij.mat (parameter file)
% [paramsPure{:}] = unifacSetUp(compArray)
% [gamma]=unifac(x, T, paramsPure{:})
% The path must be configured so that Matlab can find unifacAij.mat, unifacSetup.m, unifac.m.
% unifacAij.mat is loaded from within unifacSetup.m.

% add path for the routine unifac.m shared by all unifac routines.
addpath('../');

% provide the names for labelling using a cell array.
names = { 'Isopropanol' 'Water' 'Propionic Acid'}

% compArray is a cell array to specify the groups (not necessarily the 
% smallest groups for a molecule) When functional groups
% are near ends, they are usually kept with 
% the ends. See textbook for examples. compArray is the 
% listing of names and 'nu' values. Each column of 
% 'nu' values correspond to a molecule. Molecules should be entered
% in the same order as the names above. To use different groups, add or
% delete rows. To change number of compounds, add or delete columns.

% use "load(unifacAij.mat)" to load the group table into the workspace.
% Then open the 'groupList' table in the array editor to browse.

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
x = [0.3 0.6 0.1];
% normalize to assure sum to one.
x = x/sum(x)
TCelsius = 23

% The unifac function can be called whenever x or T changes.
% The pure component cell array is passed along with x and T.
% x should be a one row vector. A corresponding one row vector of gammas
% is returned.
[gamma] = unifac(x,TCelsius, paramsPure{:})

% ver 1.02 improved documentation regarding unifacAij.mat, unifac.m path.
% ver 1.01 improved documentation.