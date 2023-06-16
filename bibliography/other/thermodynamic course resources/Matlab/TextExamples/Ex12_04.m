function Ex12_04
% this is Example 12.4
% use of van Laar and Scatchard Hildebrand to fit methanol + benzene

% 1 - methanol
% 2 - benzene
% data - by column, sets of x(1), y(1), T(C) 
data = [
0	0	80.1;
0.026	0.267	70.67;
0.05	0.371	66.44;
0.088	0.457	62.87;
0.164	0.526	60.2;
0.333	0.559	58.64;
0.549	0.595	58.02;
0.699	0.633	58.1;
0.782	0.665	58.47;
0.898	0.76	59.9;
0.973	0.907	62.7;
1	1	64.7];
% 79 methanol
% 42 benzene
propsRow = [79 42]; % rows from ../Props/props.mat for components
antoineRow = [2  32]; %rows from ../Psat/antoineTable.mat

% ****** extract data for ease of manipulation 
xd = data(:,1);
% create matrix of x values, with pairs in each row
xd = [xd 1-xd];
y1d = data(:,2);
TKd = data(:,3) + 273.15;

ndata = size(data,1); %number of data rows

% load the properties
db = load ('../Props/props.mat'); props = db.props; clear db;
% set path for van Laar
addpath('../gammaModels/');

%solubility parameter is in the column 10
del = [props{propsRow(1),10} props{propsRow(2),10}]; 

%densities in column 7, MW in column 8
v = [props{propsRow(1),8}/props{propsRow(1),7} props{propsRow(2),8}/props{propsRow(2),7}];

%****** fit azeotrope
x = [0.614 0.386]; %azeotrope point to fit
P = 760; % mmHg
T = 58.3 + 273.15; %K
vmix = x*v'; %mixture volume
phi = x.*v/vmix;
R = 8.3145;
Psat = [591.3 368.7];

disp('Fitting azeotrope');
disp('Scatchard-Hildebrand');
[kij fval exitflag] = fzero(@calck, [-1 1]) % Scatchard-Hildebrand
disp('van Laar')
gamma = P./Psat;
lng = log(gamma); % next fit van Laar
AvL = [lng(1)*(1+x(2)*lng(2)/x(1)/lng(1))^2 lng(2)*(1+x(1)*lng(1)/x(2)/lng(2))^2]


function obj = calck(kij)
   gamma = exp((v.*(1-phi).^2)*((del(1)-del(2))^2+2*kij*del(1)*del(2))/R/T);
   obj = (x.*gamma)*Psat' - 760;
end

%calculate resulting T-x-y diagram
x = 0:0.05:1;
n = length(x);
x = [x' 1-x'];
addpath '../Psat';
[name A B C] = AntoineGet(antoineRow); %use function in ../Psat

for i=1:n
    [TK(i) fval exitflag] = fzero(@(T)bubT(T,x(i,:),kij), [300 360]);
    y(i)=y1;
end
    function obj = bubT(T,xo,kij) %objective for Scatchard-Hildebrand
        vmix = xo*v';
        phi = xo.*v./[vmix vmix];
        gamma = exp((v.*(1-phi).^2)*((del(1)-del(2))^2+2*kij*del(1)*del(2))/R/T);
        Psat = 10.^(A-B./(T-273.15 + C));
        Pcalc = (xo.*gamma)*Psat';
        y1 = xo(1)*gamma(1)*Psat(1)/Pcalc;
        obj = Pcalc - 760;
    end

plot(xd(:,1),TKd,'o'); %plot expt x vs T
hold on
plot(y1d, TKd,'^'); % plot expt y vs T
curve(1) = plot(x(:,1),TK,'-.'); % plot calc x vs T for Scatchard
plot(y,TK,'-.'); % plot calc y vs T for Scatchard


for i=1:n
    [TK(i) fval exitflag] = fzero(@(T)bubT2(T,x(i,:)), [300 360]);
    y(i)=y1;
end
    function obj = bubT2(T,xo,kij) %objective for van Laar
        [gamma(1) gamma(2)] = vanLaar(AvL(1),AvL(2), xo(1), xo(2));
        Psat = 10.^(A-B./(T-273.15 + C));
        Pcalc = (xo.*gamma)*Psat';
        y1 = xo(1)*gamma(1)*Psat(1)/Pcalc;
        obj = Pcalc - 760;
    end

curve(2) = plot(x(:,1),TK,'-'); % plot calc x vs T for van Laar
plot(y,TK,'-'); % plot calc y vs T for van Laar

xlabel('x,y methanol')
ylabel('T(K)')
legend(curve,'Scatchard Hildebrand','van Laar')

end

% ver 1.02 2/7/13 update for new props.mat
% ver 1.01 6/5/12 use struct var for 'load'
