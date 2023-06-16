% example 17.7 multiple reactions, example at one T.
T = 473;
Ka1 = exp(-96865/8.314*(1/T - 1/473.15)+3.8205);
Ka2 = exp(24155/8.314*(1/T - 1/298.15)+6.9156);
obj = @(xi)[4*xi(1)^3-Ka1*(1-xi(1)-2*xi(2))*(1+2*xi(1))^2;
           xi(2)^2 - Ka2*(1-xi(1)-2*xi(2))^2];
disp('Example at 473K')       
xi = fsolve(obj,[1;0])

% build output table like on pg 662
disp('-----'); disp('Build table like pg 662')
outputLabel = {'T(K)'; 'Ka1'; 'Ka2';'xi(1)'; 'xi(2)'; 'y1'; 'y2'; 'y3'; 'y4'; 'y5'};
output=[];
for T=473:20:573
   Ka1 = exp(-96865/8.314*(1/T - 1/473.15)+3.8205);
   Ka2 = exp(24155/8.314*(1/T - 1/298.15)+6.9156);
   obj = @(xi)[4*xi(1)^3-Ka1*(1-xi(1)-2*xi(2))*(1+2*xi(1))^2;
       xi(2)^2 - Ka2*(1-xi(1)-2*xi(2))^2];
   [xi fval exitflag] = fsolve(obj,[1;0], optimset('Display','off'));
   if (exitflag < 1) fprintf('Convergence problem at T = %g exitflag = %d\n',T,exitflag); end
   y1 = (1-xi(1)-2*xi(2))/(1+2*xi(1)); y2= xi(1)/(1+2*xi(1));
   y3 = 2*xi(1)/(1+2*xi(1)); y4=xi(2)/(1+2*xi(1)); y5 = y4;
   output= [output [T;Ka1;Ka2;xi(1);xi(2);y1;y2;y3;y4;y5]]; % build matrix of numerical results
end
out = [outputLabel num2cell(output)] % convert 'output' to cell array before cat