function RaoultFL
% vectors for Antoine Coefficients
% all values set in code, none from workspace.

%-------------input section
% Constants
names = {'pentane', 'heptane'};

% Tc(K), Pc(MPa), w in rows for each component
PureProps = [
	469.7 3.369 0.249; % pentane
	540.3 2.736 0.349; %heptane
	];
	
T = 60+273.15 % Kelvin
P = 0.1 % MPa
z = [0.5 0.5]; %overall composition

%-------------end of input section

%Extract Tc, Pc, w to make code below easier to read.
Tc = PureProps(:,1)'; %convert first column to a row vector
Pc = PureProps(:,2)'; %convert second column to a row vector
w = PureProps(:,3)'; %convert third column to a row vector
   
%normalize in case error in sum(z) = 1
z = z/sum(z);

% VOF stands for 'Vapor Over Feed', V/F.
VOFguess = 0.5;

options=optimset('Display','iter');   % Option to display output
[VOF,fval,status] = fsolve(@calcObj,VOFguess,options);  % Call optimizer

% output results
sprintf('Flash Results. T(K) = %g, P = %g',T,P)
disp(sprintf('1= %s, 2= %s',char(names(1)),char(names(2))))
format compact
VOF
z
x
y
format loose

    % ---------------------------------------
    % nested function so share variables with parent function
    function [obj] = calcObj(VOF)
    % programmed using the shortcut vapor pressure equation    
    
    Tr = T./Tc;
    Psat = Pc.*10.^(7/3*(1+w).*(1-1./Tr));
    K = Psat/P;
    x = z./(1 + VOF.*(K-1));
    y = K.*x;
    obj = sum(x) - sum(y);

    end
    %----------------------------------------
    
end %function RaoultFL

% 3/12/14 convert to V/F to be consistent with textbook



