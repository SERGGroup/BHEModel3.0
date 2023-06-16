function fitPsat
% demo of curve fitting
set(0,'format','shortg'); set(0,'formatspacing','compact');
%pairs of T(C), Psat(mmHg) for isobutanol
data = [
-9	    1;
11.6	5;
21.7	10;
32.4	20;
44.1	40;
51.7	60;
61.5	100;
75.9	200;
91.4	400;
108	    760;
];
%extract data
T = data(:,1); Pexpt = log10(data(:,2));
%initial guesses will be in array param
A = 8.5; B = 2000; C = 273;
param = [A B C];
disp('A   B   C')
[param, FVAL, EXITFLAG] = fminsearch(@calcobj,param)
%FVAL and EXITFLAG can be useful if troubleshooting is needed
disp('T  log10(Pexpt) log10(Pcalc)');
[T Pexpt Pcalc]
disp('log10(Psat) = A - B/(T + C)')
param
plot(1./(T+273.15),Pexpt,'bo');
hold on
plot(1./(T+273.15),Pcalc,'r-'); xlabel('1/T (K^{-1})'); ylabel('log_{10}P(mmHg)');

    function err = calcobj(param)
        A = param(1); B = param(2); C = param(3); %use scalars below
        Pcalc =  A - B./(T + C);
        err = sum((Pcalc - Pexpt).^2); % need scalar error for fminsearch
    end %calcobj
end