function [xR,xE,EOF] = LLEflash(z,TK)
% This function can be used to calculate LLE phase splitting
% Components 1 and 2 must be the key immiscible components.
% Put the less dense component as component 1.
% This is the key component in the extract.
% TK is in K.
% Written by Mike Dittmer and Carl Lira.

%create path to gamma models
addpath(genpath('../gammaModels'));


% Examples for UNIFAC and UNIQUAC are preprogrammed.
% option = 1 UNIFAC
% option = 2 UNIQUAC
option = 1;

switch option
    case 1,
        disp('Using UNIFAC')
        TCelsius = TK - 273.15;
        %   UNIFAC - set structure in compArray
        %       example MEK - water
        compArray = {
            'CH3CO'     1       0;
            'CH2'       1       0;
            'CH3'       1       0;
            'H2O'       0       1;
            };
        [paramsPure] = unifacSetUpLLE(compArray);
        %       then use [gamma] = unifac(x,TCelsius, paramsPure{:})
        %        Results at 298.15K       MEK   water
        %                    x(beta)     0.647     0.353
        %                    x(alpha)    0.042     0.958
    case 2;
        disp('using UNIQUAC')
        %   UNIQUAC - set parameters in uniquac.m
        %               1-Butanol ~ Water
        r = [3.4543     .92 ];
        q = [3.0520  	1.4 ];
        q1 = q;
        aij = [ 0           -82.688;
            443.56      0 ;];
        tau = exp(-aij/TK);
        %        Results at 298.15K    butanol   water
        %               x(beta)     0.50      0.50
        %               x(alpha)    0.006     0.994
    otherwise,
        disp('option is not understood')
end


%Variables
%z - Vector of feed compositions *function INPUT
%T - System Temperature in Kelvin *function INPUT
%xE - Vector of extract compositions
%xF - Vector of raffinate compositions
%EOF - E/F ratio (extract over feed)
%K - Vector of K values
%gammaE - Vector of gamma values for extract
%gammaR - Vector of gamma values for raffinate
%errorK - change in K value divided by K
%tol - Acceptable value of errorK
%count - Iteration number
%ncomp - number of components

errorK = 1;
tol = 0.00001;
count = 0;
EOF = 0.5;
ncomp = length(z);

%Start Initialization
%Initialize Raffinate and Extract Compositions
xR = zeros(1,length(z)) + 1e-50; %Zeros in the composition vectors will lead to NaN solutions
xE = xR;
xR(1:2) = [0.02 0.98];
xE(1:2) = [0.98 0.02];

%Calculate K-values
switch option
    case 1,
        gammaR = unifac(xR,TCelsius, paramsPure{:});
        gammaE = unifac(xE,TCelsius, paramsPure{:});
    case 2,
        gammaR = uniquac(xR, r, q, q1, tau);
        gammaE = uniquac(xE, r, q, q1, tau);
    otherwise,
        disp('option not understood')
end
K = gammaR./gammaE;

%End Initialization

%Start Main Loop
while errorK > tol
    
    %fzero solves the objective function(Eq 10.17) for VOF (called EOF here)
    obj = @(EOF)sum(z.*(1-K)./(1+EOF*(K-1))); % configures objective function
    EOFnew = fzero(obj,EOF); %Solves objective function, uses last EOF as initial guess
    EOF = EOFnew;

    %Calculate new compositions
    xR = z./(1+EOF*(K-1));  %Eq 10.15
    xE = xR.*K;
    %Normalize compositions
    xR = xR/sum(xR);
    xE = xE/sum(xE);

    %Calculate new K values
    switch option
        case 1,
            gammaR = unifac(xR,TCelsius, paramsPure{:});
            gammaE = unifac(xE,TCelsius, paramsPure{:});
        case 2,
            gammaR = uniquac(xR, r, q, q1, tau);
            gammaE = uniquac(xE, r, q, q1, tau);
        otherwise,
            disp('option not understood')
    end

    Knew = gammaR./gammaE;
    % echo results (use % to disable)
    fprintf('count = %d, EOF = %.4f\n', count, EOF);
    fprintf('     xR          xE       gammaR     gammaE       K\n')
    for i=1:ncomp
        fprintf('%10.4g %10.4g %10g %10g %10g %d\n', xR(i),xE(i),...
            gammaR(i), gammaE(i),  Knew(i), i);
    end
    disp('----')


    %Calculate change in K, determines if K has converged
    errorK = sum(abs(Knew - K)./(K));
    K = Knew;

    %Break out of loop if it takes too long
    count = count + 1;
    if count > 500
        error('Could not converge');
    end

    %K values will become 1 if the flash results in a one-phase product
    if ((K(1) < 1+tol) && (K(1) > 1-tol))
        error('Feed results in one phase');        
    end
end

end
%End Main Loop

% ver 1.02 converted to EOF rather than ROF. Improved iteration output.