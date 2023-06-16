function Ex17_9
% example 17.9
db=load('../Props/props.mat'); % read in as structure props
props = db.props; clear db;% convert to cell props
P = 100; %bar
Tin = 400; %K
no = [1/2 3/2 0]; % feed moles, N2, H2, NH3
nu = [-1/2 -3/2 1]; % stoich numbers
DeltaHoR = -51413; % J/mol
% look up heat capacities
% 60 N2; 59 H2; 154 NH3
id = [60 59 154]; num = length(id);
names = props(id,2);
Cp = cell2mat(props(id,15:18)); %matrix of Cp parameters
Hin = findH(Tin);

[T fval exitflag] = fsolve(@rxnadiab,600)

    function obj = rxnadiab(T)
        %equilibrium
        Ka = 0.0417659*exp(51413/8.314*(1/T - 1/600));
        M = sqrt(27)/4*P*Ka;
        xi = 1- sqrt(1-M/(1+M));
        n = no + nu*xi;
        %energy balance, molar enthalpies
        Hout = findH(T);
        obj = no*Hin - n*Hout - xi*DeltaHoR;
    end
    function H = findH(T)
        TR = 600; %K
        H = Cp(:,1)*(T-TR) + Cp(:,2)./2*(T^2-TR^2) + ...
            Cp(:,3)./3*(T^3 - TR^3)+ Cp(:,4)./4*(T^4-TR^4);
    end
end

% ver 1.02 2/7/13 update for new props.mat
% ver 1.01 6/5/12 use struct var for 'load'

