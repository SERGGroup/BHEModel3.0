function Ex11_07
% this is Margules 2param Bubble T
% all values set in code, none from workspace.
addpath(genpath('../Psat'));
addpath(genpath('../gammaModels'));

%Constants
% ------------------input section
id = [32 3]; %Antoine coeff id in database
A12 = 1.2947;
A21 = 1.8373;
Punits = 'mmHg';
Tunits = 'C';
P = 760; %mmHg
x = [0.5 0.5]; %composition
Tguess = 50; %temperature guess

%----------------------End of input section
format compact

[gamma(1) gamma(2)] = Marg2P(A12,A21, x(1), x(2)); %Margules gamma does not depend on T

[names A B C] = AntoineGet(id);
   
options=optimset('Display','iter');   % Option to display output
[T,fval,status] = fsolve(@calcObj,Tguess, options);  % Call optimizer

disp(sprintf('T(%s) = %g', Tunits, T))
x
y
sum(y)
status
format loose

    function [obj]= calcObj(T)
    Psat = 10.^(A - B./(T + C));
    % note that gammas that depend on T must be calculated here.
    K = gamma.*Psat./P;
    y = K.*x;
    obj = 1 - sum(y);
    end
end