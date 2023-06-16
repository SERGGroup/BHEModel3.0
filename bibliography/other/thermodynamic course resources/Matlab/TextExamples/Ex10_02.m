function OBJEB = Ex10_02(T)
% Adiabatic flash
% Guess exit temp T (Celsius) in workspace, Energy balance objective function
% is returned.
% A reasonable answer has 0<L/F<1.
% The OBJEB increases when T increases.

%-------------input section
% Constants
TF = 60; % Feed temp in C
Tref = 298.15; % Ref T in K
addpath ../Psat;
% methanol - 1
% ethanol - 2
[names A B C] = AntoineGet([2 3]);

%Boiling temperature in K
Tb = [64.6 78.5] + 273.15; %K
%Heat of vaporization
Hvap = [35.21 38.58] * 1000; %J/mol
% component heat capacity coefficients in rows.
% Enter 0 for fourth constant if it does not exist.
% Coefficients for Cp in J/molK.
CpL = [111.7	-4.264E-01	1.090E-03       0.;
       281.6	-1.435	2.903E-03       0.];
CpV = [21.15	0.07092	2.587E-05	-2.852E-08;
        9.014	0.2141	-8.390E-05	1.373E-09];

P = 40; % mmHg
z = [0.5 0.5]; %overall composition

%-------------end of input section

% convert Temperature to Kelvin
TFK = TF + 273.15; % Feed Temp in K
TK = T + 273.15; % Exit Temp in K

%normalize in case error in sum(z) = 1
z = z/sum(z);

LOFguess = 0.5;

options=optimset('Display','iter');   % Option to display output
[LOF,fval,status] = fsolve(@calcObj,LOFguess,options);  % Call optimizer
disp(sprintf('LOF = %g, OBJ = %g, status = %g', LOF, fval, status))

% output results
disp(sprintf('******\nFlash Results. FeedT(C) = %g, T(C) = %g, P(mmHg) = %g',TF,T,P))
disp(sprintf('1= %s, 2= %s',char(names(1)),char(names(2))))
format compact
LOF
z
x
y
K = y./x

% calculate enthalpies of the streams
HLi = [ Hcalc(CpL(1,:), Tref, TK) Hcalc(CpL(2,:),Tref, TK)];
HL = x*HLi';
HV1i = [ Hcalc(CpL(1,:), Tref, Tb(1)) Hcalc(CpL(2,:), Tref, Tb(2))]; 
HV2i = [ Hcalc(CpV(1,:), Tb(1), TK) Hcalc(CpV(2,:), Tb(2), TK)];
HV = y*[ HV1i + Hvap + HV2i ]';
HFi = [ Hcalc(CpL(1,:), Tref, TFK) Hcalc(CpL(2,:), Tref, TFK)];
HF = x*HFi';

%Check Energy balance.
OBJEB = ((1-LOF)*HV + LOF*HL - HF)/1000.;

format loose

    % ---------------------------------------
    % nested function so share variables with parent function
    function [obj] = calcObj(LOF)
    Psat = 10.^(A - B./(T + C));
    K = Psat/P;
    x = z./(K + LOF*(1-K));
    y = K.*x;
    obj = sum(x) - sum(y);

    end
    %----------------------------------------

function [H] = Hcalc(Cp,T1,T2)
% this function calculates enthalpy change given T1, T2  and the Cp
% constants. The units of Cp determin the units of H.
% The form of the polynomial as given by Elliott and Lira.
Cpc = Cp(1) + Cp(2)*T2 + Cp(3)*T2^2 + Cp(4)*T2^3; %included for convience of spot checking Cp.
H = Cp(1)*(T2 - T1) + Cp(2)/2*(T2^2-T1^2) + Cp(3)/3*(T2^3-T1^3) + Cp(4)/4*(T2^4-T1^4);
end % function Hcalc

end %function RaoultFL 



