function Rxnadiat
% example for adiabatic formation of ammonia
% limited by equilibrium
% 1/2 N2 + 3/2 H2  =  NH3
nu = [-0.5 -1.5  1];

% at 298.15K
DG298 = -16.4013; %kJ/mol
DH298 = -45.940; %kJ/mol
Tref = 298.15; %K
Tguess = 700.7;  % K
P = 100; % bar
no = [0.5 1.5  0]; %initial moles
Tin = 298.15; %K, inlet T

CpN2 = [31.2,-0.0136,2.68e-05,-1.17e-08];
CpH2 = [27.1,0.00927,-1.380e-05,7.65e-09];
CpNH3 = [27.3,0.0238,1.71e-05,-1.19e-08];
CpA = [CpN2(1) CpH2(1) CpNH3(1)];
CpB = [CpN2(2) CpH2(2) CpNH3(2)];
CpC = [CpN2(3) CpH2(3) CpNH3(3)];
CpD = [CpN2(4) CpH2(4) CpNH3(4)];

R = 0.00831447; %kJ/mol.K

% calculate inlet molar enthalpies
Hin = [enth(CpN2, Tin, Tref) ...
    enth(CpH2, Tin, Tref) ...
    enth(CpNH3, Tin, Tref)];
%convert to extensive flow
Hin = no*Hin';

DCpA = nu*CpA';
DCpB = nu*CpB';
DCpC = nu*CpC';
DCpD = nu*CpD';

[T fval exit] = fzero(@ebal, Tguess);
exit
T

function obj = ebal(T)
     %calculate K
     J=DH298+(-DCpA*298.15-DCpB/2*298.15^2-DCpC/3*298.15^3-DCpD/4*298.15^4)/1000;
     I=(DG298/298.15-J/298.15+(DCpA*log(298.15)+DCpB/2*298.15+DCpC/6*298.15^2+DCpD/12*298.15^3)/1000)/R;
     lnK=-(J/T+(-DCpA*log(T)-DCpB/2*T-DCpC/6*T^2-DCpD/12*T^3)/1000)/R-I;
     K = exp(lnK);
     M = K*P*1.29903;
     %calculate reaction coordinate
     xi = 1-sqrt(1-M/(1+M));
     %calculate moles of each specie
     n = [1.5*(1-xi) .5*(1-xi) xi];
     % calculate inlet molar enthalpies
     Hout = [enth(CpN2, T, Tref) ...
            enth(CpH2, T, Tref) ...
            enth(CpNH3, T, Tref)];
     %convert to extensive flow
     Hout = n*Hout';
     obj = Hin - Hout -xi*DH298*1000;
     disp(sprintf('T = %f, Ka = %g, obj = %f',T,K,obj))
end
    function H = enth(Cp, T, Tref)
        H = Cp(1)*(T-Tref) + Cp(2)/2*(T^2 - Tref^2) ...
            + Cp(3)/3*(T^3 - Tref^3) + Cp(4)/4*(T^4 - Tref^4);
    end
end        

