% Script file to calculate the fugacity of a mixture based on the
% Peng-Robinson Equation of State.
% This script mimics PrFug.xlsx but can handle any number of components

% In order to operate this script, the file must be located in a folder
% adjacent to the "Props" folder in the progpack directory.  The 'addpath'
% command in the code assumes this is where the file is located.  The
% command is going up one folder level in the directory, then entering the
% "Props" folder, and accessing the property database in this folder.
% Thus, if this file is not located in a folder adjacent to the "Props"
% folder, an error will occur.  File placement is crucial for proper script
% operation.  If errors continually occur, check file placement.

% Logan Matthews and Carl T. Lira
% CHE 321 Honors Option, Logan Matthews, Spring 2012

%*******************************************************
%The following lines are data to be inputted by the user

%Use propsTableBrowse (in props folder) to view the props table
%Select the table row number of your components, not the component id!
%Use a row vector only.
clear;
propsRows = [8, 11, 12];

%Input Mole Fraction Compositions. Use a row vector only in the same
%order as propsRows.
zi = [.8168 .1501  .0331];

%Input kij information as if moving across the rows of a an upper triangular matrix,  
%table, omiting the diagonals. For example, example for five components,
%the elements are entered in the order
%kij = [k12 k13 k14 k15 k23 k24 k25 k34 k35 k45]. The matrix below shows
%how the vector indices relate to row and column numbers.
%     0     1     2     3     4
%     0     0     5     6     7
%     0     0     0     8     9
%     0     0     0     0    10
%     0     0     0     0     0
% For three components, kij = [k12 k13 k23]
kij = [7.6e-4 1.71e-3 ... %make vector look like a matrix for readability
              6.1e-4];
                                  
%Input temperature in K
T = 323.15; %K
%Input Pressure in MPa
P = .2; %MPa
%End of input parameters
%*******************************************************

%Basic error checking to ensure user input was done correctly
zi = zi/sum(zi); %normalize in case user makes an error. Force sum(zi) to 1.
nComp = length(propsRows);
fprintf('\n***************\n%d components\n', nComp)

%check length of kij
if (length(kij) ~= nComp*(nComp-1)/2)
    fprintf('Binary interaction vector is size %d, expecting %d',length(kij),nComp*(nComp-1)/2)
    return;
end
if (length(zi) ~= nComp)
    fprintf('zi is size %d, expecting size %d.  Please double check input data. \n \n', length(zi),nComp);
    return;
end

addpath(genpath('../Props')); %add path to database -- Do not change unless
%file is moved from original position
if(~exist('props','var')) % check if props database is loaded.
    db = load ('props.mat'); %load as structured variable db.
    props = db.props; clear db; %transfer to cell array and clear db.
end

%Read data from the props table
for i = 1:nComp
    %Error check to ensure proper withdrawal of compound i data
    err = 0;
    if (isempty(props{propsRows(i),2}))
        fprintf('Compound %d name not found.\n',i)
        err = 1;
    else
        %Read the name of compound i
        namesA{i} = props{propsRows(i),2};
    end
    if(isempty(props{propsRows(i),4}))
        fprintf('%s Critical Temperature not found.\n',namesA{i})
        err=1;
    else
        %Read the critical temperature of compound i
        Tc(i) = props{propsRows(i),4};
    end
    if(isempty(props{propsRows(i),5}))
        fprintf('%s Critical Pressure not found.\n', namesA{i})
        err=1;
    else
        %Read the critical pressure of compound i
        Pc(i) = props{propsRows(i),5};
    end
    if(isempty(props{propsRows(i),6}))
        fprintf('%s Acentric factor not found.\n', namesA{i})
        err=1;
    else
        %Read the eccentric factor of compound i
        w(i) = props{propsRows(i),6};
    end
    if(err == 1)
        disp('Terminating. Check props database row and folder path.')
        return;
    end
end

%Calculate reduced temperatures and pressures
Tr = T./Tc;
Pr = P./Pc;
R = 8.314472; %MPa.cm^3/mol.K

%Calculate parameters for the Peng-Robinson EOS for pure fluids
b = 0.0777960739*R.*Tc./Pc; %Equation 7.16
ac = 0.4572355289.*(R.*Tc).^2./Pc; %Equation 7.16
kappa = 0.37464+1.54226.*w-0.26992.*w.^2; %Equation 7.17
alpha = (1+kappa.*(1-sqrt(Tr))).^2; %Equation 7.17

%Calculate dimensionless parameters
A = ac.*alpha.*P/(R.*T).^2; %Equation 7.21
B = b.*P./R./T; %Equation 7.22

%Construct matrix representing A.
%Amat(i,j) = Amat(j,i) is the parameter for pair i,j.
Amat = zeros(nComp, nComp);
kijmat = zeros(nComp,nComp);
%the maxtrix elemenents are stored as a packed lower triangular matrix.
%Though the description above was given as to read across the matrix,
%the order for storing an upper triangular is not as easy to understand 
%in this introductory application. The order for storing the kij gives the
%transpose of the lower triangular packed storage.
for j = 1:nComp %row
    for i = j:nComp % column, i>=j
        % determine index for kij in packed vector
        if i==j
            kijmat(i,j) = 0; % diagonals are not stored and always 0
        else
            kijmat(i,j) = kij(i + (2*nComp - j)*(j-1)/2 - j);
            kijmat(j,i) = kijmat(i,j);
        end
        Amat(i,j) = sqrt(A(i)*A(j))*(1-kijmat(i,j));
        Amat(j,i) = Amat(i,j);
    end
end

%Calculate the sum of the matrix above to satisfy the sum on page 585
AMix = zi * Amat * zi';
%Calculate the B of mixing using formula on page 585
BMix = sum(zi .* B);

%Calculate coefficients for the cubic equation to solve
a2 = -(1 - BMix);
a1 = AMix-3*BMix^2-2*BMix;
a0 = -(BMix*AMix-BMix^2-BMix^3);
%Store the coefficients in an array
polyCo = [1 a2 a1 a0];
%Use the roots solver to solve the cubic equation for Z
Ztemp = roots(polyCo);
Z = [];
%determine indices for roots are real.
index = find(imag(Ztemp)== 0);
%collect the real values of Z.
Z = real(Ztemp(index));
Z = sort(Z, 'descend');
rootRegion = length(Z);

%Calculate the volume of the mixture using Z
V = Z * R * T / P;
%Clear the array to store fugacity
fug = [];
%find all fugacities for the mixture
fug = P * exp((Z - 1) - log(Z-BMix) - AMix /BMix/2.8284*log((Z+2.4142*BMix)./(Z-.4142*BMix)));

%Clear two matrices that will store the component fugacity coefficents and
%component fugacities. Each column is for a component, and each
%row is for a root.
compFugCo = [];
compFug = [];
%loop through each root
for i = 1:rootRegion
    %calculate the component coefficients
    compFugCoI = exp(B/BMix*(Z(i)-1)-ones(1,nComp)*log(Z(i)-BMix)-AMix/BMix/2.8284*log((Z(i)+2.4142*BMix)/(Z(i)-.4142*BMix))*((2*zi*Amat)/AMix-B/BMix));
    %Use the coefficients to calculate the component fugacities
    compFugI = P .* zi .* compFugCoI;
    %Add the current calculations to the matrices storing all of the
    %component fugacities and fugacity coefficients for the current root
    compFugCo = [compFugCo; compFugCoI];
    compFug = [compFug; compFugI];
end

%Output results with proper formatting
fprintf('     Tc       Pc    w    Name\n')
for i = 1:nComp
    fprintf('%8g %8g %5g  %-12s\n',Tc(i),Pc(i),w(i),namesA{i})
end
fprintf('kij matrix\n')
format shortg; format compact
kijmat
%uncomment next statements to view matrix and/or vector
%Amat
%B

fprintf('\nT = %5g K, P = %5g MPa\n', T, P)
fprintf('zi = [ ')
for i = 1:nComp
    fprintf('%5g ',zi(i))
end
fprintf(']\n')
%Output root region
fprintf('Mixture is in %d root region\n', rootRegion);

%Output properties by rows for each component
%Each column is a different root
switch rootRegion
    case 1
        formatSpec = '%12g %-s\n';
    case 2
        formatSpec = '%12g %12g %-s\n';
    otherwise
        formatSpec = '%12g %12g %12g %-s\n';
end
fprintf(formatSpec, Z','Z')
fprintf(formatSpec, V', 'V(cc/mol)')
fprintf(formatSpec, fug','f_mix(MPa)')
fprintf('component fugacity coeff\n')
for i=1:nComp
    fprintf(formatSpec,compFugCo(:,i),namesA{i})
end
fprintf('component fugacity (MPa)\n')
for i=1:nComp
    fprintf(formatSpec,compFug(:,i),namesA{i})
end

% ver 1.0 6/5/12 initial release
