% Example 14_03
% water + MEK immiscibility by DeltaG mixing curve
%create path to gamma models
addpath(genpath('../gammaModels'));

names = {'water','methylethylketone'};
compArray = {
   'CH3CO'     0       1;
   'CH2'       0       1;
   'CH3'       0       1;
   'H2O'       1       0;
   };

T = 298; % K
TCelsius = T - 273.15;

[paramsPure] = unifacSetUpLLE(compArray);
R = 8.314; % J/molK

% (using a very small number instead of zero avoids blowup of log)
x = [1e-10 .01:.01:0.99 1-1e-10];
x = [x' 1-x']; %matrix of x1, x2 in columns
size = length(x);

% set up array for Gtotal
DGmixoRT = zeros(1,size);

% loop over composition
for i = 1: size
   gamma = unifac(x(i,:),TCelsius, paramsPure{:});
   DGmixoRT(i)= x(i,:)*log((x(i,:).*gamma)');
end

plot (x(:,1),DGmixoRT)
xlabel('x_{water}')
ylabel('\DeltaG_{mix}/RT')

    