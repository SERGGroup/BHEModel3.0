%   Peng Robinson Property Calculator

%This script solves the Peng Robion cubic using T(K), P(MPa) from workspace.
%This script illustrates
%   how a script leaves all variables in the workspace
%   solving polynomials for roots
%   sorting vectors
%   eliminating imaginary roots

% Developed by Carl Lira for use with
% Introductory Chemical Engineering Thermodynamics, by
% J. Richard Elliott and Carl T. Lira, Prentice-Hall, 2nd ed., 2013.
% See http://chethermo.net for links to 
% download an original. To be redistributed only in un-modified form, 
% and only for educational use. The program developers have no liability 
% for the use of the program or results. Copyright 2010-13, Carl T. Lira.

%*****************Input section**************************   
    %Set the ID from ../Props/props.mat. Use 'run ../Props/propsTableBrowse' to preview.
    propsRow = 67; %enter the row number from the table, not the ID number
%*****************End of input section********************	
	
%{
    Specify T(K) and P(MPa) in the workspace by entering into 
    the workspace numerical values calculated from the critical 
    properties (Tc, Pc) of the compound you have specified.
    To try out the calculator, try the pairs
        T = 0.7 * Tc, P = 0.5 * Pc
        T = 0.7 * Tc, P = 0.2 * Pc
    and T = 1.5 * Tc, P = 0.5 * Pc
    (where Tc and Pc values are given below)
    For propane these pairs of values are:
        T = 259; P = 2.13; 
        T = 259; P = 0.85;
        T = 555; P = 2.13;

    Then run the script by entering >> 'run Preos'.

    Note how a script writes all the values to the workspace.

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

if(err == 1)
   disp('Terminating. Check props database row and folder path.')
   return;
end

% set output to eliminate spaces between lines
format compact

sprintf('%s  Tc(K)= %5.1f Pc(MPa)= %5.3f w = %5.3f',name,Tc, Pc, w)

if (not(exist('P')) | not(exist('T')))
    if not(exist('P')) disp('Please specify P(MPa) and re-run preos.') 
    end
    if not(exist('T')) disp('Please specify T(K) and re-run preos.') 
    end
    return
end

% echo the values of name, T and P
sprintf('T(K)= %f P(MPa)= %f', T, P) 

R = 8.314472; %MPa.cm^3/mol.K
b = 0.0777960739*R*Tc/Pc;
ac = 0.4572355289*(R*Tc)^2/Pc;
kappa = 0.37464+1.54226*w-0.26992*w^2;

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
Zvals = roots([1 a2 a1 a0])

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
 
 % display results
 sprintf('Z= %f  %f', Z) 
 V = Z*R*T/P; %cm^3/mol
 sprintf('V(cm^3/mol)= %f  %f', V) 
  % a constant that will be useful in formulas
  sq = sqrt(2);
  
 %calculate the fugacity
 fug = P*exp(Z-1-log(Z-B)-A/B/2/sq*log((Z + (1+sq)*B)./(Z + (1-sq)*B)));
 sprintf('fugacity (MPa)= %f  %f', fug) 
 %calculate the enthalpy departure.
 Hdep=R*T*(Z-1-A/B/2/sq*log((Z+(1+sq)*B)./(Z+(1-sq)*B))*(1+kappa*sqrt(Tr)/sqrt(alpha)));
 sprintf('Hdep (J/mol)= %f  %f', Hdep) 
 
 % restore default format
 format loose
 
 % ver 1.03 2/7/13 update for new props.mat.
 % ver 1.02 6/5/12 converted 'load' to use structured variable method.
 % ver 1.01 4/9/12 modified for added short name column in props.mat
 