function [Z, H, S, U, phi, info] = PreosProps(Tref, Pref, T, P, match)
% Property calculator/solver for the Peng Robinson equation of state, 
% This file is intended to be distributed with PreosPropsMenu.m and PreosPropsMenu.fig
% Edit this file to set the fluid to be studied before using the function. 

% This file is distributed with 'fminsearch' enabled and 'lsqnonlin'
% disabled. The program is more robust with 'lsqnonlin', but some shared
% computer labs may have limited licenses for the optimization toolbox.

% Developed by Carl Lira for use with
% Introductory Chemical Engineering Thermodynamics, by
% J. Richard Elliott and Carl T. Lira, 2nd ed. Prentice-Hall, 2013.
% See http://chethermo.net for links to 
% download an original. To be redistributed only in un-modified form, 
% and only for educational use. The program developers have no liability 
% for the use of the program or results. Copyright 2009-2013, Carl T. Lira.

%-----
% The following variables(constants) are global and should never change
% values in any function after being set in the main program.
global Tc Pc w Cp R HdepRef SdepRef ac b kappa
%-----

%******************** Input section, can be edited **************
    %Set the ID from ../Props/props.mat. Use 'run ../Props/propsTableBrowse' to preview.
    propsRow = 62; %enter the row number from the table, not the ID number

%*****************End of input section********************	

%{
    To run the function from the menu, enter
	> run PreosPropsMenu
	To run the function from the command line, 
	set values for input variables Tref(K), Pref(MPa), T, P
	and 'match' (described below). Then type in the command window
	> [Z H S U phi info] = PreosProps(Tref, Pref, T, P, match);
    Or copy it from the function definition in an editor window
    and paste it into the command window and add a semicolon.
    (Be sure to include the trailing semicolon to supress the long display 
     of unformatted results after the list of formatted results).
%}

% Most variables have obvious definitions from their names.
% The info variable is a cell array that passes back several pieces of information:
% info{1}  = T
% info{2} = P
% info{3} = exitflag
% info{4} = name of compound
% info(5) = obj.

% This function
%   illustrates use of functions
%   solving and sorting polynomials
%   use of vector math and element-by-element vector math
%   use of fzero 1D solver
%   anonymous functions
%   passing input variables from the workspace to the function
%   passing values from the function back to the workspace
%   control of program options using input vector and conditional statements.

%{ 
   Peng Robinson Property Calculator for U,H,f,S
   To use this routine:   
    Set the critical properties below 
    and the heat capacity and save the m-file.    

    Specify reference temperature(K) Tref and reference pressure(MPa) Pref 
    in the workspace. 
    Convenient values are Tref = 298.15, Pref = 0.1.
    
    Specify T in K and P in MPa in the workspace.
    Suggested introductory settings are given below.

    Specify the array 'match' as described next.
        (use match = [ 0 0 0 0] to skip matching).

    The array 'match' is used to control how properties are calculated.
    When match(1) = 0, the specified T and P are used to calculate U,H,f,S
    If matching a property value is desired (match(1) ~= 0), T and/or P are 
    initial values and either T or P is adjusted to match a specifed value
    for a property. 
    Specify match=[i j k l] in the command window as follows
     i - run match?  |  j - adjust    | k - value       | l - match which root?
     ---------------------------------------------------------------------------
     0 - no matching |  (ignored for  |  (ignored for   |  (ignored for i = 0,5,6)
     1 - match U     |  i = 0)        |   i = 0,5,6     |   1 match largest Z root
     2 - match H     |  1 - adjust T  |  value to match |   2 match smallest Z root
     3 - match S     |  2 - adjust P  |                 |  
     4 - match V     |  3 - T and P   |                 |
     5 - saturation  |                |                 |       
         (match f's) |                |                 | 
     6 - user

    By performing the above steps, the workspace is initialized.

    Suggested introductory calculations
    For the specified fluid, set using the critical properties
            T = 0.7 * Tc, P = 0.2 * Pc, match = [ 0 0 0 0]; calculate props
             (for CO2 the values are T = 213; P = 1.5; match = [0 0 0 0])
            T = 0.7 * Tc, P = 0.07 * Pc, match = [ 4 2 0 0]; solve for Psat
             (for CO2 the values are T = 213; P = 0.5; match = [4 2 0 0])
             (fugacities should be matched.)
            Then try to match a value of S using the previous output
            to select a new value to specify.

    Following is a list of user-defined routines found below.
    CalcProps  - calculate the themodynamic properties.
    SortZ - determine the stable Z value using fugacity values.
    CalcPhi - calculate the fugacity values.
    FindMatch - contols program flow depending on specification of 'match'.
    PRsolveZ - solves the Peng-Robinson eq for Z values at specified T,P.
    CalcHdep - calculate the enthalpy departure for the specified Z values.

    Matlab routine fzero is used in the main routine to call FindMatch
    
%}

addpath(genpath('../Props')); %add path to database	

if(~exist('props','var')) % check if props database is loaded.
    db = load('../Props/props.mat'); %load as structured variable.
    props = db.props; clear db; %transfer to props and then clear struct.
end

name = props{propsRow,2};
Tc = props{propsRow,4};
Pc = props{propsRow,5};
w = props{propsRow,6};
Cp = [props{propsRow,15} props{propsRow,16} props{propsRow,17} props{propsRow,18} ];

err = 0;
if(isempty(name))
	disp('Compound name not found.')
	err=1;
end
if(isempty(Tc))
	disp('Critical Temperature not found.')
	err=1;
end

if(isempty(Pc))
	disp('Critical Pressure not found.')
	err=1;
end

if(isempty(w))
	disp('Acentric factor not found.')
	err=1;
end

if(numel(Cp)==0)
	disp('Heat capacity not found.')
	err=1;
end

if(err == 1)
   disp('Terminating. Check props database row and folder path.')
   return;
end

% set remaining global constants
R = 8.314472; %MPa.cm^3/mol.K
b = 0.0777960739*R*Tc/Pc;
ac = 0.4572355289*(R*Tc)^2/Pc;
kappa = 0.37464+1.54226*w-0.26992*w^2;
% end of setting constants

% calculate reference state properties
% solve for the values of Z at reference state
[Zref, Tr, A, B, alpha] = PRsolveZ(Tref,Pref);

% find the correct value of Z
[Zref, phiRef] = SortZ(Zref, A, B); 

% find the values of H, S, departure at refstate
[HdepRef] = CalcHdep(Zref, A, B, Tref, alpha);
SdepRef = HdepRef/Tref - R*log(phiRef);
% ---------------End of Ref State calculation

exitflag = 0; % initialize so that value can be passed back even if solver not used.
obj = 0;

% section for matching a state property
if match(1) ~= 0   % only optimize if variable match(1) is non-zero
    switch match(2)
       case 1
        % adjust T
        disp(sprintf('Adjusting T for fixed P(MPa)= %g',P))
        % optimset is a matlab library routine to set options.
        options=optimset('Display','iter');   % Option to display output
        % An anonymous function is used to capture exisitng state
        % variables and adjust only one of them. 
        % See Optimization Toolbox> Tutorial> Examples that use standard
        % agorithms> Avoiding global variables via anonymous...
        % fzero is a matlab library routine for non-linear solving
        [T,obj,exitflag] = fzero(@(T) FindMatch(Tref, Pref, T,P,match),T,options);  % Call optimizer
        if abs(obj) > 1E-3
            exitflag = -20; % found discontinuous point
            disp(sprintf('\nFound discontinuous point. Custom exitflag code -20.'))
        end
        if exitflag < 0
            disp(sprintf('ABNORMAL TERMINATION. Not converged.'))
        end
        disp(sprintf('\nT(K)= %g P(MPa) = %g obj= %g exitflag= %d',T,P,obj, exitflag))
        disp(sprintf('See Matlab fzero help to understand standard exitflag codes.\n'))
       case 2
        % adjust P
        disp(sprintf('Adjusting P for fixed T(K)= %g',T))
        % see comments in case 1 immediately above for optimset and fzero.
        options=optimset('Display','iter');   % Option to display output
        [P,obj,exitflag] = fzero(@(P) FindMatch(Tref, Pref, T,P,match),P,options);  % Call optimizer
        if abs(obj) > 1E-3
            exitflag = -20; % found discontinuous point
            disp(sprintf('\nFound discontinuous point. Custom exitflag code -20.'))
        end
        if exitflag < 0
            disp(sprintf('ABNORMAL TERMINATION. Not converged.'))
        end
        disp(sprintf('\nT(K)= %g P(MPa) = %g obj= %g exitflag= %d',T,P,obj,exitflag))
        disp(sprintf('See Matlab fzero help to understand exitflag codes.\n'))
        case 3
        % adjust T and P
        disp(sprintf('Adjusting T and P using initial guess T(K)= %g, P(MPa)= %g',T,P))
        % see comments in case 1 immediately above for optimset and fzero.
        options=optimset('Display','iter');   % Option to display output
        % Use lsqnonlin if you have optimization toolbox
        % Function userObj is used for user created objective.
        % The 'obj' returned by userObj must be a vector
        TnP = [T P]; %variable TnP holds 'T and P' in an array
        %[TnP, resnorm, obj, exitflag] = lsqnonlin(@(TnP)userObj(Tref, Pref, TnP,match),TnP);
        % Use fminsearch if you don't have optimization toolbox
        % The 'obj' returned by userObj must be a scalar
        [TnP, obj, exitflag] = fminsearch(@(TnP)userObj(Tref, Pref, TnP,match),TnP);
        T = TnP(1);
        P = TnP(2);
        otherwise
         error('Terminating. match(2)=%d not understood', match(2))
         return;
    end % switch match(2)
else
  % generate value of user Objective function
  TnP = [T P];
  obj = userObj(Tref, Pref, TnP, match);
end
% end of section for matching a state property

% calculate properties at the values of T, P
[Z, H, S, U, phi, Hdep, Sdep, Udep] = CalcProps(Tref, Pref, T,P); 

% display results
% note: nesting the sprint within a disp supresses the 'ans ='
format compact; %eliminate spaces between lines
% note: nesting the sprint within a disp supresses the 'ans ='
disp(sprintf('%s  Tc(K)= %g Pc(MPa)= %g w = %5.3f',name,Tc, Pc, w))
disp(sprintf('Tref= %g, Pref= %g, Zref= %g  %g', Tref, Pref, Zref)) 
disp(sprintf('T(K)= %g P(MPa)= %g', T, P))
disp(sprintf('match = %g %g %g %g', match))
disp(sprintf('Z=              %g\t  %g', Z))
V = Z*R*T/P; %cm^3/mol
disp(sprintf('(H-Hig) (J/mol)=      %g\t  %g', Hdep))
disp(sprintf('(U-Uig) (J/mol)=     %g\t  %g', Udep))
disp(sprintf('(S-Sig) (J/molK)=     %g\t  %g', Sdep))
disp(sprintf('V(cm^3/mol)=    %g\t  %g', V) )
disp(sprintf('H (J/mol)=      %g\t  %g', H))
disp(sprintf('U (J/mol)=     %g\t  %g', U))
disp(sprintf('S (J/molK)=     %g\t  %g', S))
disp(sprintf('fugacity (MPa)= %g\t  %g', phi*P))
disp(sprintf('Objective Function= %g\t %g',obj))
        if exitflag < 0
            disp('ABNORMAL TERMINATION. Not converged.See above for exitflag code.')
        end

if match(1) ~= 0
    disp(sprintf('\nManually change the value of T/P to save in workspace.'))
end
format loose; % turn on default spacing
% build info to return
info = {T P exitflag name obj};

end
% -------- End of main function

% *************************************
function obj = userObj(Tref, Pref, TnP, match)
% This function for a user objective function when match(1) = 6.
% See function used when 1 < match(1) < 5.
T = TnP(1);
P = TnP(2);
[Z, H, S, U, phi] = CalcProps(Tref,Pref, T,P);
switch match(4)
 % select the specified root
 case 1
 % index for largest root
 index = find(Z == max(Z));
 case 2
 % index for smallest root
 index = find(Z == min(Z));  
end
% Use Z(index), H(index), S(index), U(index) phi(index) in the obj
% calculation.
% Use an ARRAY obj if called from lsqnonlin when you have optimization
% toolbox
%    obj(1) = (Z(index)*R*T/P)/88.46 - 1;
%    obj(2) = U(index)/(-5330) - 1;
% Use a SCALAR obj if called from fminsearch when you do not have
% optimization toolbox.
   obj1 = 0;
   obj2 = 0;
   obj = obj1^2 + obj2^2;
end


% *************************************
function obj = FindMatch(Tref, Pref, T,P,match)
% This routine is used by fzero to generate objective function when 
% a property is to be matched; it is ignored otherwise.
% Note: this function is used when 1 < match(1) < 5.
% When match(1) = 6, then userObj is used instead of FindMatch.
global R;
[Z, H, S, U, phi] = CalcProps(Tref,Pref, T,P); 
switch match(4)
 case 1
 % index for largest root
 index = find(Z == max(Z));
 case 2
 % index for smallest root
 index = find(Z == min(Z));  
end
switch match(1)
  case 1
  obj = U(index) - match(3);
  case 2
  obj = H(index) - match(3);
  case 3
  obj = S(index) - match(3);
  case 4
  obj = Z(index)*R*T/P - match(3);
  case 5
  % index is not used
  if length(phi) < 2
      error('Fatal condition, only one root found while iterating on Psat. Terminating.')
      return;
  end %if  
  obj = phi(1)/phi(2) - 1;
  otherwise
   error('Terminating. match(1)=%d is not understood.', match(1))
  return
end
end % function FindMatch  

%****************************
function [Z, H, S, U, phi, Hdep, Sdep, Udep] = CalcProps(Tref,Pref, T,P)
% This is the main section to calculate departures and ideal gas contributions
global Cp R HdepRef SdepRef

%-Calculate departures
%-- first find roots.
[Z, Tr, A, B, alpha] = PRsolveZ(T,P);
[phi] = CalcPhi(Z, A, B);
% find the departure values of H, S
[Hdep] = CalcHdep(Z, A, B, T, alpha);
Udep = Hdep - (Z-1)*R*T;
Sdep = Hdep/T - R*log(phi);
%- End of Departure Calculations

%- Calculate ideal gas contribution
HidGas = Cp(1)*(T - Tref) + Cp(2)/2*(T^2 - Tref^2) + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4);
SidGas = Cp(1)*log(T/Tref) + Cp(2)*(T - Tref) + Cp(3)/2*(T^2 - Tref^2) + Cp(4)/3*(T^3 - Tref^3) - R*log(P/Pref);
%-End of ig contribution

%- Calculate properties
H = Hdep + HidGas - HdepRef;
S = Sdep + SidGas - SdepRef;
U = H - Z*R*T;
end

%**************************************
function [phi] = CalcPhi(Z, A, B)
   % sq is a constant that will be useful in formulas
   sq = sqrt(2);
   % perform element-by-element processing of Z vector
   phi = exp(Z-1-log(Z-B)-A/B/2/sq*log((Z + (1+sq)*B)./(Z + (1-sq)*B)));
end

%***********************************
function [Hdep] = CalcHdep(Z, A, B, T, alpha)
global Tc R kappa
 sqrt2 = sqrt(2);
 Tr = T/Tc;
 % perform element-by-element processing of Z vector
 Hdep=R*T*(Z-1-A/B/2/sqrt2*log((Z+(1+sqrt2)*B)./(Z+(1-sqrt2)*B))*(1+kappa*sqrt(Tr)/sqrt(alpha)));
end

%************************************
function [Z, phi]= SortZ(Z, A, B)
% Determine which root has the minimum fugacity, 
% using only the values > B
index = find(Z>B);
phi = CalcPhi(Z(index),A,B);

% if there is more than one root, choose the most stable
if length(index)>1
   % find the index with the minimum fugacity
   index=find(phi == min(phi));
end % if length
Z = Z(index);
phi = phi(index);
end

%*******************************************
function [Z, Tr, A, B, alpha] = PRsolveZ(T,P)
%Determine the coefficients of the cubic and solve for Z keeping
%at most two real roots.
%The global constants below must be previously set.
%global constants that should not be changed.
global Tc Pc w R ac b kappa

%calculate conditions
Tr = T/Tc;
Pr = P/Pc;
alpha = (1+kappa*(1-sqrt(Tr)))^2;
%calculate dimensionless parameters
A = ac*alpha*P/(R*T)^2;
B = b*P/R/T;

%determine cubic coefficients and solve cubic
a2 = -(1-B);
a1 = A-3*B^2-2*B;
a0 = -B*(A-B-B^2);

% This function finds the real roots.
% Zvals holds the real/imaginary results.
Zvals = roots([1 a2 a1 a0]);

%determine indices for roots are real and greater than B.
index = find(imag(Zvals)== 0);

%collect the real values of Z.
if length(index)>1
    %execute this if more than one root
    %collect the real values
    Zreal=real(Zvals(index));
    %accept the largest and smallest real values because
    %the center is always unstable when three exist
    Z = [max(Zreal), min(Zreal)];
else
    %execute this if there is just one real root
    Z = real(Zvals(index));
end
end

% ver 1.44 - 3/7/14 minor fix to insert commas in returned vectors
% ver 1.43 - 2/7/13 modified for new version of props.m.
%                   disabled lsqnonlin by default
% ver 1.42 - 6/5/12 modified use 'load' with structured variable.
% ver 1.41 - 4/9/12 modified for revised props.mat with short name column
% ver 1.4 - 3/15/12 modified to match V
% ver 1.31 - 1/19/10 added use of props database
% ver 1.3 - 10/16/09 modified to introduce userObj
% ver 1.2 - 2/17/09 added error checking for discontinuous function for fzero
% ver 1.1 - 2/22/09 added 'info' to return 'info' to Menu and workspace. Added menu.
% ver 1.0 - initial release 2008