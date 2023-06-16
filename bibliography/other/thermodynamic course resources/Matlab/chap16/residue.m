function xstore = residue
clc
% insert id rows of ../Psat/AntoineTable.mat
% Double click in Current Directory to load into workspace, 
% then double click in workspace to load into array editor.
% benzene ethanol methanol
id = [32 3 2];

% Set composition along residue curve. 
% Run multiple times ad different comps to add lines.
xstart = [0.3 0.6 0.1];
% normalize
xstart = xstart/sum(xstart);

% programmed to use UNIQUAC
% Benzene  Ethanol Methanol
r = [3.1878	2.1055	1.4311];
q = [2.400	1.972	1.432];
% the q1 is for the 'Anderson' UNIQUAC variation for increased accuracy for alcohols. 
% Anderson, T.F., Prausnitz, J.,M., Ind. Eng. Chem. Process Des. Dev. 17, 1978, 552-561.
% for standard UNIQUAC leave q1 = q;
q1 = q;
% be very careful with signs and units. This code is distributed
% using the convention tau(i,j) = exp(-a(i,j)/T) where T is in K.
% Some publications use RT instead of T. Aspen uses the opposite sign.
aij = [   0	337.7386	559.4861; 
        -42.6567	0	124.4848;
        -27.2253	-91.2264	0;];
 

if(length(id)~=3)    
    disp('The residue calculator is designed for three components.')
    return
end %if
addpath '../Psat';
addpath '../gammaModels';
[names A B C] = AntoineGet(id);
names

x = xstart;

% %dlnnL is equivalent to per cent change in liquid moles.
% %dlnL is also equal to dnL/nL
dlnnL = -0.15;

TCelsius = 70; % initial guess of bubble temperature
P = 760; % mmHg
%options=optimset('Display','iter');
% set xstore
xstore = [];
Tstore = [];


% calculate upwards
while min(x) > 0.005  
    %loop runs until practically only one of the components is gone.
    x = x/sum(x); %normalized x, in case sum of x is not 1.
    %display(x) %displays x-vector in output

    xstore = [xstore; x];

    % use fzero method to get bubble point temperature
    [TCelsius,fval] = fzero(@calcObj,TCelsius);

    Tstore = [Tstore; TCelsius];

    %TCelsius and Psat values are obtained from the objective function

    %Calculate dx vector
    %dx is the change in mole fraction corresponding to change in dlnnL
    dx = (K-1).*x*dlnnL; %term by term multiplication
    %update x vector
    x = x+dx;
end %while

% now loop upwards
dlnnL = -dlnnL;
x = xstart;
 
while min(x)>0.005 
%loop runs until practically only one of the components is gone.
x = x/sum(x); %normalized x, in case sum of x is not 1.
%display(x) %displays x-vector in output

xstore = [x; xstore];

% use fzero method to get bubble point temperature
[TCelsius,fval] = fzero(@calcObj,TCelsius);

Tstore = [TCelsius; Tstore];

%TCelsius and Psat values are obtained from the objective function

%Calculate dx vector
%dx is the change in mole fraction corresponding to change in dlnnL
dx = (K-1).*x*dlnnL; %term by term multiplication
%update x vector
x = x+dx;
end %while

%plot diagonal
plot([0 1],[1 0],'color','black');
xlabel(char(names(1)))
ylabel(char(names(2)))
hold on;

%plot data
plot(xstore(:,1),xstore(:,2))

xstore =[Tstore xstore];

%Note that 'hold' is left on so more runs will add more lines.

    function obj = calcObj(TCelsius)
        %The bubble temperature occurs when the sum of y values is 1.
        %Properties
        %Calculate vector of vapor pressures
        Psat = 10.^[A - B./(TCelsius+C)];  %mm Hg
        % matrix of tau
        TK = TCelsius + 273.15;
        tau=exp(-aij./TK);
        gamma = uniquac(x,r,q,q1,tau);
        %Calculate vector of K ratios
        K = gamma.* Psat/P;
        %Calculate vector of vapor phase composiiton y values
        y = x.*K;
        %Check if y values sum to 1.
        obj = sum(y) - 1.;
    end
end

%ver 1.01 - 6/5/12 made function name consistent with file name
