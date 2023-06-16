function [paramsPure] = unifacSetUpLLE(compArray)
% This function is to set up the UNIFAC calculation
% Before calling this function set up the compArray
% such as shown here in columns for Isopropanol, Water,
% and Propionic Acid:
%
% compArray is a table of group names from unifacAij.mat
% along with the 'nu' values for each compound in rows.
% 

db = load ('unifacAijLLE.mat');
unifacAijLLE = db.unifacAijLLE; clear db;

% determine number of data rows in table (two header rows)
% These are hard coded and must be reevaluated if unifacAij.mat
% is changed.
nHeaderRows = 2;
nDataRows = size(unifacAijLLE,1);
nComp = size(compArray,2)-1;
nGroups = size(compArray, 1);

% extract the numerical portion of compArray for easier calculations
% The 'cat' function used to dump the cell array to a number array will 
% create a row array of values that is reshaped to stay as the numerical matrix.
compNArray = reshape(cat(1,compArray{:,2:nComp+1}),nGroups,nComp);

% set up storage for molecular and group volume and surface params.
r = zeros(1,nComp);
q = zeros(1,nComp);
R = zeros(1,nGroups);
Q = zeros(1,nGroups);
aij = zeros(nGroups,nGroups);

% Now some housekeeping to make the routine user-friendly for typos
% in the group names. If the user has a typo we need to know.
% Vector 'found' will be updated to 1 for each row when name is found.
found = zeros(1,nGroups);
% keep a listing of the rows of the table that are used because they will
% be needed to build the aij matrix
rowsFound = zeros(1,nComp);

% loop through the rows of compArray
for compRow=1:nGroups
    % loop through the data table looking for the same group name
    for dataRow=3:nDataRows
        if(strcmp(compArray{compRow,1},unifacAijLLE{dataRow,2})==1)
            % the name has been found
            found(compRow)=1; % store a 'flag' to tells us this row is found
            % store the group volume and surface area
            R(compRow) = unifacAijLLE{dataRow,3};
            Q(compRow) = unifacAijLLE{dataRow,4};
            % save a record of the rows that will need to be used for
            % building aij.
            rowsFound(compRow) = dataRow;
            % increment the molecular volume and surface area
            r = r + compNArray(compRow,:)*R(compRow); 
            q = q + compNArray(compRow,:)*Q(compRow); 
            break; % break out of for loop, no reason to search further
        end %if
    end % for dataRow
end % for compRow

notfound = found - 1; %check for any missing rows, indicated by -1 entries in notfound

% for user-friendly program, give a listing of the text from compArray
% where names were not found.
if sum(notfound)< 0  % will execute if any entries are < 0
    for compRow=1:nGroups % loop to find which rows are not found.
        if notfound(compRow) < 0
            sprintf('Group %s from row %d of compArray is not in database. Check input for typo.\n', char(compArray(compRow,1)), compRow)
        end %if notfound(i)
    end % for compRow
    return; 
end %if notfound

% fill aij matrix
for i = 1:nGroups
    for j = 1:nGroups
        aij(i,j) = unifacAijLLE{rowsFound(i),rowsFound(j)+2};
    end
end

% calculate group mole fractions.
% Each compound will be in column and the groups will
% appear in rows in the same order as compArray.
% X(i,j) will be the mole fraction of the group from
% row i in compound j.

% the kron function is combined with a ones matrix
% to build a matrix using copies of a vector. Here the 
% assembled matrix is used for element-by element division.
Xpure = compNArray ./ kron(sum(compNArray,1),ones(nGroups,1));

% sum(X*Q) is needed for surface area fractions
% Because Q is a row it looks 'backwards'.
sumXQ = Q*Xpure;

% Group surface area fractions in a matrix.
% Subscripts work as in Xpure.
THETApure = Xpure.*kron(Q',ones(1,nComp))./kron(sumXQ,ones(nGroups,1));

% pack into a cell array for convenience in passing to other routines.
paramsPure = {r, q, R, Q, nComp, nGroups, compNArray, Xpure, THETApure, aij};

end


% ver 1.02 4/20/13 improved documentation.
% ver 1.01 6/5/12 modified to use 'load' with stuctured variable.
