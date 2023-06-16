function [paramsPure] = unifacSetUp(compArray)
% This function is to set up the UNIFAC calculation
% Before calling this function set up the compArray
% such as shown here in columns for Isopropanol, Water,
% and Propionic Acid:
%
% compArray is a table of group names from unifacAij.mat
% along with the 'nu' values for each compound in rows.
% 
% compArray = { 
%     'CH3'       2       0       1;
%     'CH2'       0       0       1;
%     'CH'        1       0       0;
%     'OH'        1       0       0;
%     'H2O'       0       1       0;
%     'COOH'      0       0       1;
%     };

db = load ('unifacAij.mat');
groupList = db.groupList; unifacAij = db.unifacAij; clear db;

% Determine number of header and data rows in cell Arrays unifacAij and groupList
% Header rows are hard coded and must be reevaluated if unifacAij.mat
% or groupList are changed.
nHeaderRowsG = 1;
nHeaderRowsA = 2;
nDataRowsA = size(unifacAij,1);
nDataRowsG = size(groupList,1);
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
% Also, we need to track if the main group is found in the Aij table.
foundAij = zeros(1,nGroups);
rowAij = zeros(1,nGroups);

% Allocate storage for a listing of the main groups
mainGroup = zeros(1,nComp);

% loop through the rows of compArray
for compRow=1:nGroups
    % loop through the subgroup data table looking for the same subgroup group name
    for gRow=(nHeaderRowsG+1):nDataRowsG
        if(strcmp(compArray{compRow,1},groupList{gRow,2})==1)
            % the name has been found
            found(compRow)=1; % store a 'flag' to tells us this row is found
            % store the group volume and surface area
            R(compRow) = groupList{gRow,5};
            Q(compRow) = groupList{gRow,6};
            % save a record of the main group number
            mainGroup(compRow) = groupList{gRow,3};
            % find the main group location in the unifacAij table
            for aRow=(nHeaderRowsA+1):nDataRowsA
                % For extracting the numerical value, the element in the 
                % cell array must be converted to a
                % number. This is handled by addressing the cell array
                % element.
                temp = unifacAij{aRow, 1};            
                if ( mainGroup(compRow) == temp)
                    foundAij(compRow)=1;
                    rowAij(compRow)=aRow;
                    break;
                end 
            end %aRow
            % increment the molecular volume and surface area
            r = r + compNArray(compRow,:)*R(compRow); 
            q = q + compNArray(compRow,:)*Q(compRow); 
            break; % break out of for loop, no reason to search further
        end %if
    end % for dataRow
end % for compRow

notfound = found - 1; %check for any missing rows, indicated by -1 entries in notfound
notfoundAij = foundAij - 1; %check for any missing Aij rows;

% for user-friendly program, give a listing of the text from compArray
% where names were not found.
if (sum(notfound)< 0 || sum(notfoundAij) < 0) % will execute if any entries are < 0
    for compRow=1:nGroups % loop to find which rows are not found.
        if notfound(compRow) < 0
            sprintf('Subgroup %s from row %d of compArray is not in database. Check input for typo.\n', char(compArray(compRow,1)), compRow)
        end %if notfound(i)
        if notfoundAij(compRow) < 0
            sprintf('Subgroup %s, Maingroup %d from row %d of compArray is not found in the Aij matrix. Check input for typo.\n', char(compArray(compRow,1)), mainGroup(compRow), compRow)
        end %if notfoundAij(i)
    end % for compRow
    return; 
end %if notfound

% fill aij matrix
for i = 1:nGroups
    for j = 1:nGroups
        aij(i,j) = unifacAij{rowAij(i),rowAij(j)};
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

% pack parameters into a cell array for convenience in passing to other functions.
paramsPure = {r, q, R, Q, nComp, nGroups, compNArray, Xpure, THETApure, aij};

end

% ver 1.03 4/20/13 improved documentation
% ver 1.02 6/5/12 modified to use struct var when using 'load'
% 6/10/09 Modified to use separate listing of Groups from the Aij matrix










        


