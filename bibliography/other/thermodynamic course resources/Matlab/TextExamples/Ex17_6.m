% example 17.6 multiple reactions
set(0,'format','shortg');set(0,'formatspacing','compact')
Ka1 = 2.667; Ka2 = 3.2;
obj = @(xi)[Ka1*(2-xi(1)-xi(2))*(1-xi(1))-(xi(1)-xi(2))*xi(1);
           Ka2*(2-xi(1)-xi(2))*(xi(1)-xi(2))-4*xi(2)^2];
disp('Results for Example 17.6')       
xi = fsolve(obj,[0;0]); xi

 
