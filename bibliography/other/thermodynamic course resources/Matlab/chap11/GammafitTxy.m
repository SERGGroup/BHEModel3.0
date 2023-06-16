function GammafitTxy

% Input definitions;
% x1 is an array with values in the range of 0<x1<1
% A12 & A21 are adjustable parameters 
% Texp,Tcalc is the temperature
% x1 is the composition of the liquid composition of component 1

%-------------------------------------------------------------------------
% Provide experimental data to match

%Data are in columns, x1, y1, T, 
Data = [
    0	0	100;
    0.0115	0.163	95.17;
    0.016	0.2115	93.4;
    0.0365	0.3655	88.05;
    0.057	0.4565	84.57;
    0.1	0.5015	82.7;
    0.1215	0.512	82.32;
    0.1665	0.5215	81.99;
    0.1895	0.5375	81.58;
    0.1935	0.532	81.75;
    0.245	0.539	81.62;
    0.2835	0.553	81.23;
    0.2975	0.554	81.29;
    0.298	0.551	81.28;
    0.3835	0.57	80.9;
    0.446	0.592	80.67;
    0.5145	0.6075	80.38;
    0.559	0.6255	80.31;
    0.646	0.6645	80.15;
    0.6605	0.6715	80.16;
    0.6955	0.6915	80.11;
    0.765	0.737	80.23;
    0.809	0.7745	80.37;
    0.8725	0.834	80.7;
    0.9535	0.9325	81.48;
    1	1	82.25];

Aij = [0 0]; % intial parameter guesses.
id = [6 44]; %ids for isopropanol and water Antoine constants
Pexpt = 760; % mmHg

R = 8.31447; %J/mol-K
x1 = Data(:,1);
ndata = length(x1);
y1expt = Data(:,2); 
Texp = Data(:,3);

%-------------------------------------------------------------------------
addpath(genpath('../Psat'));
addpath(genpath('../gammaModels'));

% Component 1
[names A B C] = AntoineGet(id)

x2 = 1-x1;              % calculate x2
x  = [x1  x2];         % create a 2 column array of x1 & x2

% disable the opimzation call to see the effect of A12 and A21 by
% trial/error.
% Use lsqnonlin if you have optimization toolbox
% The 'obj' returned by calcError must be a vector
Aij = lsqnonlin(@calcError,Aij)
% Use fminsearch if you don't have optimization toolbox

A12 = Aij(1);
A21 = Aij(2);

% generate Txy using fitted parameters.
x1c = 0:0.05:1 + 1E-50;
x2c = 1-x1c + 1E-50;
ncalc = length(x1c);
ycalc=zeros(ncalc);
for j = 1:ncalc
    Tcalc(j) =fzero(@Tbub,40);
    [Gamma1Calc Gamma2Calc] = Marg2P(Aij(1), Aij(2), x1c(j), x2c(j));
    Psat = 10.^(A - (B./(Tcalc(j)+C)));
    Pcalc = x1c(j)*Gamma1Calc*Psat(1) + x2c(j)*Gamma2Calc*Psat(2);
    y1c(j)= x1c(j)*Gamma1Calc*Psat(1)/Pcalc;
end

% plot results
plot(x1,Texp,'o');
hold on;
plot(x1c,Tcalc);
plot(y1c,Tcalc);
xlabel('x-y'); ylabel('T(C)')

    function obj = Tbub(T)
        % This function is used to generate Txy after fitting
        [Gamma1Calc Gamma2Calc] = Marg2P(Aij(1), Aij(2), x1c(j), x2c(j));
        Psat = 10.^(A - (B./(T+C)));
        Pcalc = x1c(j)*Gamma1Calc*Psat(1) + x2c(j)*Gamma2Calc*Psat(2);
        obj = Pcalc-Pexpt;
    end

%-------------------------------------------------------------------------
% Define objective function for parameter fitting
% minimize sum(Pexpt - Pcalc)^2

    function obj = calcError(Aij)
        % This function is used to fit parameters to experimental data.
        %system compositions [x1 x2]
        A12 = Aij(1);
        A21 = Aij(2);
        obj = [];
        for i = 1:ndata
        [Gamma1Calc Gamma2Calc] = Marg2P(Aij(1), Aij(2), x1(i), x2(i));
        Psat = 10.^(A - (B./(Texp(i)+C)));
        Pcalc = x1(i)*Gamma1Calc*Psat(1) + x2(i)*Gamma2Calc*Psat(2);
        % use this obj for lsqnonlin if you have optimization toolbox
        obj(i) = Pcalc - Pexpt;
        % use this obj for fminsearch
        % obj = obj + (Pcalc - Pexpt).^2;
        end       
        % uncomment next line to echo progress to screen
        % [Aij obj]
    end
end

% Revision summary
% 1.00 Initial release 12/15/12
