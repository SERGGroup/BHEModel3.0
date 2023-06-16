function Kcalc
% This function calculates the equilibrium constant.

% The function can be modified to use fzero.

% Developed by Carl Lira, Copyright 2010-2012.

% First column of id is the props row from ../Props/props.mat.
% Second column is the stoichiometric number.

%{ 
example CO + 2H2 = MeOH
   row      61   59      79
    nu      -1   -2      1

id = [ 61  -1;
       59 -2;
       79  1];
%}

id = [ 61  -1;
       59 -2;
       79  1];
   
db = load ('../Props/props.mat');
props = db.props; clear db;

disp('*****Ka calculator******')

R = 0.008314472; %kJ/mol.K
n = size(id,1); %number of species
useCp = 1; % 1-yes, 0-no, will be set to zero if shortcut van't Hoff is to be used.

disp('nu     DHf298     DGf298   name')

for i=1:n
    CpA(i)=0; CpB(i)=0; CpC(i)=0; CpD(i)=0;
    nu(i) = id(i,2); %variable copy for clarity of programming
    name(i) = props(id(i,1),2);
    if isempty(props{id(i,1),13})
        disp(sprintf('DHf not found for %s. Terminating', name{i}))
        return
    else
        DHf(i) = props{id(i,1),13};     
    end
    if isempty(props{id(i,1),14})
        disp(sprintf('DGf not found for %s. Terminating', name{i}))
        return
    else
        DGf(i) = props{id(i,1),14};       
    end
    if isempty(props{id(i,1),15})
        useCp = 0; %reset flag 
    else
        CpA(i) = props{id(i,1),15};
    end

    if ~isempty(props{id(i,1),16})
        CpB(i) = props{id(i,1),16};
    end
    if ~isempty(props{id(i,1),17})
        CpC(i) = props{id(i,1),17};
    end
    if ~isempty(props{id(i,1),18})
        CpD(i) = props{id(i,1),18};
    end
    disp(sprintf('%3d %10g %10g %s',nu(i), DHf(i), DGf(i),name{i}))
end
disp('       CpA        CpB        CpC        CpD   name')
for i=1:n
    disp(sprintf('%10g %10g %10g %10g %s',CpA(i), CpB(i), CpC(i), CpD(i),name{i}))
end

DH298 = nu*DHf';
DG298 = nu*DGf';
disp('DH298(kJ/mol)  DG298(kJ/mol)')
disp(sprintf('%10g %10g',DH298, DG298))

msg = 'Would you like to use the short-cut van''t Hoff?';
if useCp == 0 %confirm to use short-cut if Cp coefficients missing
    msg = strcat('Heat capacity coefficients missing for a component. ',msg);
    useCp = menu(msg,'Yes','No');
    useCp = useCp-1;
    if useCp ==1
        disp('Unable to continue. Cp constants missing.')
        return
    end
end

if useCp == 1 % give choice to use shortcut if Cp all valid
    useCp = menu(msg, 'Yes', 'No');
    useCp = useCp -1;
end

DCpA = 0; %set deltaCp values
DCpB = 0;
DCpC = 0;
DCpD = 0;

if useCp == 1
    DCpA = nu*CpA';
    DCpB = nu*CpB';
    DCpC = nu*CpC';
    DCpD = nu*CpD';
    fprintf('Using full van''t Hoff\n')
else
    fprintf('Using shortcut van''t Hoff\n')
end

J=DH298+(-DCpA*298.15-DCpB/2*298.15^2-DCpC/3*298.15^3-DCpD/4*298.15^4)/1000;
I=(DG298/298.15-J/298.15+(DCpA*log(298.15)+DCpB/2*298.15+DCpC/6*298.15^2+DCpD/12*298.15^3)/1000)/R;

T = 298.15; % default value to start
while T > 0
    DHT = J+(DCpA*T+DCpB/2*T^2+DCpC/3*T^3+DCpD/4*T^4)/1000;
    lnK=-(J/T+(-DCpA*log(T)-DCpB/2*T-DCpC/6*T^2-DCpD/12*T^3)/1000)/R-I;
    DGT = -lnK*R*T;
    K = exp(lnK);
    disp(sprintf('T(K)= %g, K= %g',T,K))
    disp(sprintf('DHT(kJ/mol)= %g, DGT(kJ/mol)= %g',DHT,DGT))
    T = input('Enter T in Kelvin. Use a negative value to quit. T = ');
end

% ver 1.04 2/7/13 modified for new version of props.mat
% ver 1.03 5/8/12 fixed handling of empty cell array properties.
% ver 1.02 4/13/12 loading props into structured variable.
% ver 1.01 4/9/12 modified to use revised props.mat with added short name
%    column. fixed an issue with id(i) instead of id(i,1).