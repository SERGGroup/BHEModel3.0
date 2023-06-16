function examples
% This function shows how to call the activity coefficient models.
% See readme.txt for more details.
set(0,'format','shortg'); set(0,'formatspacing','compact');

%*****Binary Models that can accept arrays of binary compositions ****
x1 = 0:0.1:1;
x2 = 1 - x1;
A12 = 1.;
A21 = 1.5;
[gamma1, gamma2] = Marg1P(A12,x1, x2);
fprintf('\nOne-Parameter Margules, A12 = %g\n',A12)
disp('x1 gamma1 gamma1')
[x1' gamma1' gamma2']

[gamma1, gamma2] = Marg2P(A12,A21, x1, x2);
fprintf('\nTwo-Parameter Margules, {A12,A21} = {%g, %g}\n',A12, A21)
disp('x1 gamma1 gamma2')
[x1' gamma1' gamma2']

[gamma1, gamma2] = vanLaar(A12,A21, x1, x2);
fprintf('\nvan Laar, {A12,A21} = {%g, %g}\n',A12, A21)
disp('x1 gamma1 gamma2')
[x1' gamma1' gamma2']

%******* models that accept mulicomponent compositons calculate
% only one compositon at a time.

TK = 298.15; %Kelvin
x = [0.1 0.2 0.3 0.4]; x = x/sum(x);
%binary parameters for NRTL.
%see below, tau=aij+bij/TK.

% For the parameter version in the textbook, set aij matrix to zeros.
% Diagonals are zero.
aij=[0	0	0	-1.9763;
0	0	1.817306	0.806535;
0	-4.41293	0	-2.34561;
3.3293	0.514285	3.853826	0;	
];

% Diagonals of bij are zero.
bij=[0	-252.4821	-235.2789	609.8886;
225.4756	0	-421.289	-266.533;
515.8212	1614.287	0	1290.464;
-723.8881	444.8857	-4.42868	0;
];

% symmetric matrix of non-randomness parameters.
alpha=[0	0.3	0.3	0.3;
0.3	0	0.1	0.4;
0.3	0.1	 0	0.364313;
0.3	0.4	 0.364313	0;
];

%See Aspen reference for NRTL model
%tau is truncated after second term
%alpha is truncated after first term
tau=aij+bij/TK;
gamma=nrtl(x, tau, alpha);
disp(' '); disp('NRTL')
disp('x gamma')
[x; gamma]

% **** UNIQUAC *******************
% 1 - ethanol
% 2 - tolueme
x = [0.3 0.7]; x = x/sum(x);
r = [2.1055 3.9228];
q = [1.972 2.968];
% the q1 is for the 'Anderson' UNIQUAC variation for increased accuracy for alcohols. 
% Anderson, T.F., Prausnitz, J.,M., Ind. Eng. Chem. Process Des. Dev. 17, 1978, 552-561.
% for standard UNIQUAC leave q1 = q;
q1 = q;

% be very careful with signs and units. This code is distributed
% using the convention tau(i,j) = exp(-a(i,j)/T) where T is in K.
% Some publications use RT instead of T. 
% Aspen uses the opposite sign compared to most literature, and 
% the aij in literature is equivalent to -bij in Aspen, when Apsen aij = 0.
% aij has zeros on the diagonal.
aij = [   0   -76.1573;
         438.005  0];
% matrix of tau
tau=exp(-aij./TK);
[gamma] = uniquac(x, r, q, q1, tau);
disp(' '); disp('UNIQUAC; ethanol + toluene')
disp('x gamma')
[x; gamma]

% water, methyl ethyl ketone (MEK)
r = [.92 3.2479];
q = [1.4 2.876];
q1 = q;
aij = [ 0    -2.0882;
      345.53      0;
      ];
% matrix of tau
tau=exp(-aij./TK);
[gamma] = uniquac(x, r, q, q1, tau);
disp(' '); disp('UNIQUAC; water + MEK')
disp('x gamma')
[x; gamma]

% See UnifacVLE/UnifacCaller.m
% and UnifacLLE/UnifacCallerLLE.m for examples of how to call UNIFAC.

end

% ver 1.02 4/20/13 improved documentation.
% ver 1.01 used correct set statements for shortg and compact.