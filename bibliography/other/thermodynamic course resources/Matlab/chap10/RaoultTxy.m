function RaoultTxy
% Raoult's Law Txy using Antoine Psat
% vectors for Antoine Coefficients
% all values set in code, none from workspace.

% The routine uses the same methods as RaoultBT, and the background
% of the calculation methods and variables are largely described in
% that function. 

% This routine uses
%   element-by-element vector math
%   1-D solver
%   plotting multiple data sets and labelling chart

% Developed by Carl Lira for use with
% Introductory Chemical Engineering Thermodynamics, by
% J. Richard Elliott and Carl T. Lira, Prentice-Hall, 2nd ed., 2013.
% See http://chethermo.net for links to 
% download an original. To be redistributed only in un-modified form, 
% and only for educational use. The program developers have no liability 
% for the use of the program or results. Copyright 2008-2013, Carl T. Lira.
addpath(genpath('../Psat')); %add path to Antoine constants

%**************Beginning of input section, editable
% id values from rows of ../Psat/AntoineTable.mat
id = [27 31];
%***************End of input section******************

[names A B C] = AntoineGet(id);
names

% Specify units for Antoine Equation constants.
Tunits = 'C';
Punits = 'mmHg';

%variables
Tguess = 50; %initial guess
% set P basis for bubble T calculation
P = 1000; % same units as Antoine's equation

% Create composition vectors that span the composition range.
x1 = 0:0.1:1;
x2 = 1 - x1;

% Create rows of composition pairs in X by tranposing the 
% composition vectors.
X = [x1' x2'];
% count number of rows
Nrows = length(x1);

% find the bubble pressure for each row of X (each pair of compositions)
for i = 1:Nrows

    options=optimset('Display','iter');   % Option to display output
    [T,fval,status] = fzero(@calcObj,Tguess);  % Call optimizer

    % store the results for later plotting
    Tstore(i) = T;
    % y values are available because a nested function is used.
    Ystore(i,:) = y;

end %for X loop

%print table of results
%'sprintf' is nested in 'disp' to create tighter output format
disp(sprintf('1 = %s, 2 = %s',char(names(1)),char(names(2))))
disp(sprintf('P = %g %s', P, Punits))
disp(sprintf('x1\ty1\tT(%s)',Tunits))
for i = 1:Nrows
    disp(sprintf('%g\t%g\t%g',X(i,1),Ystore(i,1),Tstore(i)))
end

%plot results
% plot x vs. T with a default solid line
plot(X(:,1),Tstore)
% set plot limits for xaxis
xlim([0 1])
hold on % add to plot
% plot y vs. T on same plot with dashed line
plot(Ystore(:,1),Tstore,'--')
% add chart labels
xlabel(['x - y ' char(names(1))])
ylabel(['T(' Tunits ')'])
title([char(names(1)) ' + ' char(names(2)) ' at P= ' num2str(P) ' ' Punits])
hold off
% end of main function

% ------------------------------------------------------------
%beginning of nested objective function(shares vars with main function)
    function obj = calcObj(TCelsius)

        %properties
        % vector of vapor pressures.
        % Use element-by-element vector math
        Psat = 10.^(A - B./(TCelsius+C));  %mm Hg
        % vector of K ratios
        K = Psat/P;
        % use the active row of X to calculate the corresponding y's
        y = X(i,:).*K;
        % objective function 
        obj = sum(y) - 1.;
    end
% ------------------------------------------------------------

end

% ver 1.1 Added use of ../Psat/AntoineTable 1/19/10 CTL
