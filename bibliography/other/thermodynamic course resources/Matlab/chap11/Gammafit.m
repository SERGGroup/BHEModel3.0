function Gammafit

% Input definitions;
% x1 is an array with values in the range of 0<x1<1
% A12 & A21 are adjustable parameters used in matching the 
% T is the temperature
% x1 is the composition of the liquid composition of component 1

% This code is distributed with a 'pause()' statement that stops at each iteration.
% Documentation near the 'pause()' statement describes the behavior.

%-------------------------------------------------------------------------
% Provide experimental data to match

%Data are in columns, x1, y1, P, 
Data = [
         0         0   32.1000;
    0.0015    0.0254   33.8000;
    0.0111    0.1374   37.1000;
    0.0231    0.2603   42.3000;
    0.0357    0.3577   47.2000;
    0.0649    0.4604   55.0000;
    0.1168    0.5316   60.3000;
    0.1970    0.5547   62.9000;
    0.2271    0.5611   63.5000;
    0.3120    0.5659   64.4000;
    0.3958    0.5907   65.1000;
    0.4477    0.5890   65.8000;
    0.5009    0.6098   66.6000;
    0.6369    0.6462   66.9000;
    0.7542    0.7296   66.8000;
    0.8245    0.7752   65.7000;
    0.9363    0.8892   63.2000;
    1.0000    1.0000   60.7000;
];

x1 = Data(:,1);
y1expt = Data(:,2); 
Pexpt = Data(:,3);

id = [6 44]; %ids for isopropanol and water Antoine constants
T = 80; % in degrees C

%-------------------------------------------------------------------------
addpath(genpath('../Psat'));
addpath(genpath('../gammaModels'));

% Component 1
[names A B C] = AntoineGet(id)
Psat = 10.^(A - (B./(T+C)));
Psat(1) = 60.7; % override calculated values
Psat(2) = 32.1;

x2 = 1-x1;              % calculate x2
x  = [x1  x2];         % create a 2 column array of x1 & x2

% initial guess for A12 and A21
Aij = [0 0];

% disable the opimzation call to see the effect of A12 and A21 by
% trial/error.
% Use lsqnonlin if you have optimization toolbox
% The 'obj' returned by calcError must be a vector
Aij = lsqnonlin(@calcError,Aij)
% Use fminsearch if you don't have optimization toolbox
% The 'obj' returned by calcError must be a scalar
% Aij = fminsearch(@calcError,Aij)

A12 = Aij(1);
A21 = Aij(2);

% Calculate Gamma for defined parameters of A1 and A2
[Gamma1Calc Gamma2Calc] = Marg2P(A12, A21, x1, x2);

Pcalc = (x1.*Gamma1Calc)*(Psat(1)) + (x2.*Gamma2Calc)*(Psat(2));
y1 = x1.*Gamma1Calc*Psat(1)./Pcalc;

plot(x1,Pcalc)
hold on
plot(y1,Pcalc)
plot(x1,Pexpt,'s')
plot(y1expt,Pexpt,'^')
hold off

%-------------------------------------------------------------------------
% Define objective function
% minimize sum(Pexpt - Pcalc)^2

    function obj = calcError(Aij)

        %Define system compositions [x1 x2]

        A12 = Aij(1);
        A21 = Aij(2);
        % Calculate Gamma for defined parameters of A1 and A2
        [Gamma1Calc Gamma2Calc] = Marg2P(A12, A21, x1, x2);

        Pcalc = (x1.*Gamma1Calc)*(Psat(1)) + (x2.*Gamma2Calc)*(Psat(2));

        % use the obj for lsqnonlin if you have optimization toolbox
        obj = Pcalc - Pexpt;
        % use this obj for fminsearch
        % obj = sum( (Pcalc - Pexpt).^2 );
        
        % echo progress to screen
        Aij
        [x1 Pexpt Pcalc Pcalc-Pexpt]
        y1 = x1.*Gamma1Calc*Psat(1)./Pcalc;
        plot(x1,Pcalc)
hold on
plot(y1,Pcalc)
plot(x1,Pexpt,'s')
plot(y1expt,Pexpt,'^')
hold off
% The pause() statement stops execution each time it is encountered until the spacebar is pressed.
% The lsqnonlin routine evaluates the current iteration, then makes an infinitesimal step in each
% variable sucessively to evaluate the derivative of the objective function with respect to each variable.
% Then the next step is taken. The diagram will typically look identical for 3 times the pause is encountered
% when there are two parameters - (evaluate, test param1, test param2), then step. Insert a value (seconds)
% in pause to stop for a fixed amount of time - pause(2) will pause for 2 seconds.
pause()

    end
end

% Revision summary
% 1.15 3/27/13 Documentation of the 'pause' statement was added.
% 1.1 12/15/12 'A' was used for Margules parameters and Antoine
% coefficients. The duplication did not cause any numerical errors because
% the Antoine equation calculation was overridden. The Margules parameters
% were changed to Aij.