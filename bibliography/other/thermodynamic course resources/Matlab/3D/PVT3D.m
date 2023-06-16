function PVT3D
%Tom Fielitz Honors Option 1 CHE 321 Spring 2011

% the subcritical region will be generated in
% three parts: vapor, liquid, two-phase.
% All three parts will be plotted separately.

% All three parts will use common temperatures
% so that they look continuous.

% variable ranges
%figure
lowt = 5; %deg C
% The XSteam routine does not return the 
% psat correctly when exact Tc is input.
% Therefore, use a T that is close.
% critt = 374.16; %deg C
critt = 373.9; % deg C, approx critt
lowp = .01; %bar
critp = 221.2; %bar
hight = 500; % deg C
highp = 500; %bar


% number of increments for each of the three regions
numTsub = 20; %spans from lowt to critt
numPVap = 20; %spans from lowp to critp
numPLiq = 20;  %spans from lowp to highp

% declare arrays for each of the regions
% Tsub - array of T's to be used for subcritical properties
%        The same T's will be used for V, VL, L reqions.
% PVap - array of pressures to be used for vapor region
% PLiq - array of pressures to be used for liquid region
% The sub critical surfaces will be created at the grid of points
% created by combining Tsub with the P arrays.
% For the two phase region, only Tsub is needed to find Psat.
Tsub = linspace(lowt,critt,numTsub);
PVap = linspace(lowp,critp,numPVap);
PLiq = linspace(lowp,highp,numPLiq);
% set up zero valued arrays to copy temperatures from Tsub
TarrayVap = zeros(numTsub,numPVap);
TarrayLiq = zeros(numTsub,numPLiq);
TarrayLV = zeros(11,numTsub);
% set up zero valued arrays to store the pressures
ParrayVap = zeros(numTsub,numPVap);
ParrayLiq = zeros(numTsub,numPLiq);
ParrayLV = zeros(11,numTsub);
% set up zero valued arrays to store enthalpies
vV = zeros(numTsub,numPVap);
vL = zeros(numTsub,numPLiq);
vLV = zeros(11,numTsub);

% % create the subcritical vapor region
for column1 = 1:numTsub
    % create a copy of Tsub and insert it into TarrayVap
    t = Tsub(column1);
    TarrayVap(:,column1) = t;
    % determine Psat at each T
    psat = XSteam('psat_T',t);
    for row1 = 1:numPVap
        % for each of the stored pressures
        p = PVap(row1);
        if p > psat
            % to make plot display correctly, set
            % all points for this criteria 
            % to Psat and VsatV. It will
            % 'smash' them all together, but the
            % plot looks correct.
            ParrayVap(row1,column1) = psat;
            vV(row1,column1) = XSteam('vV_p',psat);
        else
            % This is the region that will create the
            % surface that will be visible.
            % Use the P and vapor V
            ParrayVap(row1,column1) = p;
            vV(row1,column1) = XSteam('v_pt',p,t);
        end %if
    end %for p
end %for t
% draw the surface
surf(vV, TarrayVap, ParrayVap);
hold all;
%
%
% To understand what has been done to this point,
% insert a breakpoint here.
% 
% create the liquid region
for column2 = 1:numTsub
    t = Tsub(column2);
    TarrayLiq(:,column2) = t;
    psat = XSteam('psat_T',t);
    for row2 = 1:numPLiq
        p = PLiq(row2);
        if p < psat
            % to make plot display correctly, set
            % all points for this criteria 
            % to Psat and HsatL. It will
            % 'smash' them all together, but the
            % plot looks correct. 
           ParrayLiq(row2,column2) = psat;
           vL(row2,column2) = XSteam('vL_p',psat);
        else
            % use the P and liquid H
            ParrayLiq(row2,column2) = p;
            vL(row2,column2) = XSteam('v_pt',p,t);
        end %if
    end %for p
end %for t
surf(vL,TarrayLiq,ParrayLiq);
%
%
% To understand what has been done to this point,
% insert a breakpoint here.
%
% create the saturated two phase region
for column3 = 1:numTsub
    t=Tsub(column3);
    TarrayLV(:,column3) = t;
    ParrayLV(:,column3) = XSteam('psat_T',t);
    vSatL = XSteam('vL_t',t);
    vSatV = XSteam('vV_t',t);
    % use 10 intervals to increment quality by 0.1
    for row3 = 0:10
        q = row3/10;   % quality
        vLV(row3+1,column3) = q * vSatV  + (1-q) * vSatL ;
    end
end
surf(vLV,TarrayLV,ParrayLV);

%
% To understand what has been done to this point,
% insert a breakpoint here.
%
%create supercritical region
%p>critp
%T<tcrit


% declare arrays

PLiqsup = [ParrayVap(:,20)' ParrayLiq(10:20,20)'];

number = length(PLiqsup);
%'number' is necessary dimension of arrays to fit desired pressure range
numTsup = number; 
numberT = round(numTsub*((hight - critt))/(critt-lowt)); 
%Sizes increments of T axis to match subcritical region

%will only create 'numberT' number of Temperatures
        %All other temperatures will be default hight to stack ontop one
        %another
numPLiqsup = number;  %spans from critp to highp
%PLiqsup = linspace(critp,highp,numPLiqsup);
TarrayVapsup = zeros(numTsup,number);
TarrayLiqsup = zeros(numTsup,number);
%set default temperature array to highest temperature so that any
%additional temperature values stack at highest temp.
for column5 = 1:number
    TarrayVapsup(:,column5) = hight;
end

Tsup = [linspace(critt,hight,numberT) TarrayVapsup(1,numberT:31)];
    %defines Tsup as the vector of superheated temperatures with equal
    %spacing over the range, then calls a row of '500's from TarrayVapsup
    %to stack the rest of the temperatures ontop of one another at the edge
    %of the graph

%Create arrays
ParrayVapsup = zeros(numTsup,number);
ParrayLiqsup = zeros(numTsup,number);
vVsup = zeros(numTsup,number);
vLsup = zeros(numTsup,number);


for column4 = 1:number
    t = Tsup(column4);
    TarrayLiqsup(:,column4) = t;
    for row4 = 1:numPLiqsup
        p = PLiqsup(row4);
        ParrayLiqsup(row4,:) = p;
        vLsup(row4,column4) = XSteam('v_pt',p,t);
    end %for p
end %for t

%add volume at critical point from sub-critical vapor region
vLsup(20,1) = vLV(1,20);
surf(vLsup,TarrayLiqsup,ParrayLiqsup);

surf(vVsup,TarrayVapsup,ParrayVapsup);
% adjust axes
axis([0 500 0 hight .01 highp]);
axis manual;
% make V axis be logarithmic
 set(gca,'xscale','log');    


%load a colormap to make shading nicer
%  db = load('CmapLira.mat')
%  colormap(db.mycmap);
%enable (if desired) once colormap is loaded into appropriate folder

caxis([-100 550]);
%Shifts color spectrum down and up so that the deep blue is below the P axis
%and the dark red is above the top of the P

% label the axis
title('Pressure-Volume-Temperature Chart for H_2O');
xlabel('Log Volume (m^3/kg)');
ylabel('Temperature (deg C)');
zlabel('Pressure (bar)');

% adjust the cammera angle for a nice appearance
cameraangle = [11973,-5832.8,3603.5];
campos(cameraangle);
hold off
end

% ver 1.01 6/5/12 converted 'load' to use structured variable method.
