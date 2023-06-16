% This script is an example that calls UNIFAC

% The unifac calculation requires the following
% files (also showing the calling method):
% unifacAijLLE.mat (parameter file)
% [paramsPure{:}] = unifacSetUpLLE(compArray)
% [gamma]=unifac(x, T, paramsPure{:})
% The path must be configured so that Matlab can find unifacAijLLE.mat, unifacSetupLLE.m, unifac.m.
% unifacAijLLE.mat is loaded from within unifacSetupLLE.m.

% add path to unifac.m shared with all unifac codes.
addpath('../');
% provide the names for labelling
names = { 'Isopropanol' 'Water' 'Propionic Acid'}

% compArray is a cell array to specify the groups (not necessarily the 
% smallest groups for a molecule) When functional groups
% are near ends, they are usually kept with 
% the ends. See textbook for examples. compArray is the 
% listing of names and 'nu' values. Each column of 
% 'nu' values correspond to a molecule. Molecules should be entered
% in the same order as the names above. To use different groups, add or
% delete rows. To change number of compounds, add or delete columns.

% use "load('unifacAijLLE.mat')" to load the database into the workspace.
% Open the table in the array editor to browse available groups.

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
[paramsPure] = unifacSetUpLLE(compArray);

% set concentration and temperature.
x = [0.8 0.1 0.1];
% normalize to assure sum to one.
x = x/sum(x)
TCelsius = 23

% The unifac function can be called whenever x or T changes.
% The pure component cell array is passed along with x and T.
% x should be a one row vector. A corresponding one row vector of gammas
% is returned.
[gamma] = unifac(x,TCelsius, paramsPure{:})

% ver 1.02 improved documentation regarding unifacAijLLE.mat, unifacSetupLLE.m path.
% ver 1.01 improved documentation