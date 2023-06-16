function pht3D
% the subcritical region will be generated in
% three parts: vapor, liquid, two-phase.
% All three parts will be plotted separately.
% All three parts will use common temperatures
% so that they look continuous.

% Created by Eric Selden and Carl Lira
% For distribution with Introductory Chemical Engineering Thermodynamics
% J.R. Elliott, C.T. Lira, 2nd ed., Prentice-Hall.

% This code requires that XSteam.m and CmapLira.mat be in the Matlab path.
% XSteam is from www.x-eng.com. CmapLira.mat is a color map to make the surface less dark.

% variable ranges
lowt = 5; %K
% The XSteam routine does not return the 
% psat correctly when exact Tc is input.
% Therefore, use a T that is close.
% critt = 374.16; %deg C
critt = 373.9; % deg C, approx critt
lowp = .01; %bar
critp = 221.2; %bar
hight = 800; %bar
highp = 1000; %bar


% number of increments for each of the three regions
numTsub = 10; %spans from lowt to critt
numPVap = 10; %spans from lowp to critp
numPLiq = 10;  %spans from lowp to highp

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
hV = zeros(numTsub,numPVap);
hL = zeros(numTsub,numPLiq);
hLV = zeros(11,numTsub);

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
            % to Psat and HsatV. It will
            % 'smash' them all together, but the
            % plot looks correct.
            ParrayVap(row1,column1) = psat;
            hV(row1,column1) = XSteam('hV_p',psat);
        else
            % This is the region that will create the
            % surface that will be visible.
            % Use the P and vapor H
            ParrayVap(row1,column1) = p;
            hV(row1,column1) = XSteam('h_pt',p,t);
        end %if
    end %for p
end %for t
% draw the surface
surf(hV, TarrayVap, ParrayVap);
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
           hL(row2,column2) = XSteam('hL_p',psat);
        else
            % use the P and liquid H
            ParrayLiq(row2,column2) = p;
            hL(row2,column2) = XSteam('h_pt',p,t);
        end %if
    end %for p
end %for t
surf(hL,TarrayLiq,ParrayLiq);
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
    hSatL = XSteam('hL_t',t);
    hSatV = XSteam('hV_t',t);
    % use 10 intervals to increment quality by 0.1
    for row3 = 0:10
        q = row3/10;   % quality
        hLV(row3+1,column3) = q * hSatV  + (1-q) * hSatL ;
    end
end
surf(hLV,TarrayLV,ParrayLV);

%
% To understand what has been done to this point,
% insert a breakpoint here.
%

%create supercritical region
%p>critp
%T<tcrit

numTsup = 10; %spans from critt to hight
numPVapsup = 10; %spans from lowp to critp
numPLiqsup = 10;  %spans from critp to highp

% declare arrays
Tsup = linspace(critt,hight,numTsup);
PVapsup = linspace(lowp,critp,numPVapsup);
PLiqsup = linspace(critp,highp,numPLiqsup);
TarrayVapsup = zeros(numTsup,numPVapsup);
TarrayLiqsup = zeros(numTsup,numPLiqsup);
ParrayVapsup = zeros(numTsup,numPVapsup);
ParrayLiqsup = zeros(numTsup,numPLiqsup);
hVsup = zeros(numTsup,numPVapsup);
hLsup = zeros(numTsup,numPLiqsup);


for column4 = 1:numTsup
    t = Tsup(column4);
    TarrayLiqsup(:,column4) = t;
    for row4 = 1:numPLiqsup
        p = PLiqsup(row4);
        ParrayLiqsup(row4,:) = p;
        hLsup(row4,column4) = XSteam('h_pt',p,t);
    end %for p
end %for t
surf(hLsup,TarrayLiqsup,ParrayLiqsup);


for column5 = 1:numTsup
    t = Tsup(column5);
    TarrayVapsup(:,column5) = t;
    for row5 = 1:numPVapsup
        p = PVapsup(row5);
        ParrayVapsup(row5,:) = p;
        hVsup(row5,column5) = XSteam('h_pt',p,t);
    end %for p
end %for t

surf(hVsup,TarrayVapsup,ParrayVapsup);
% adjust axes
axis([0 4500 0 800 .01 1000]);
axis manual;
% make P axis be logarithmic
set(gca,'zscale','log');    


%load a colormap to make shading nicer
db = load('CmapLira.mat');
colormap(db.mycmap);
caxis([-200 1000]);

% label the axis
xlabel('Enthalpy (kJ/kg)');
ylabel('Temperature (deg C)');
zlabel('Log Pressure (bar)');

% adjust the cammera angle for a nice appearance
cameraangle = [11973,-5832.8,3603.5];
campos(cameraangle);

% ver 1.01 6/5/12 converted 'load' to use structured variable method.



