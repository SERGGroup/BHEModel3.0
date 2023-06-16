% Example 17.3 reworked showing that Gibbs energy is minimized
% at equilibrium conversion. From Example 17.3, the equilibrium
% conversion is 0.784 at 1 bar and 900K.

names = {'butane','butadiene','hydrogen', 'water'};
% butane  = butadiene + H2
% steam is 'inert'
% at reaction temperature from textbook example
% Note: it is necessary to recalculate DGform if T changes.
DGform = [232.854 243.474 0 -198.204];
T = 900; % K

R = 8.314; % J/molK
P = 1; %bar
DGoRT = 1000*DGform/R/T; %convert to J/mol

% reaction coordinate array
% (using a very small number instead of zero avoids blowup of log)
xi = [1e-10 .05:.05:0.95 1-1e-10];
size = length(xi);

% set up array for Gtotal
Gtot = zeros(1,size);

% specify stoichometry
n.i = [ 1 0 0 10];
nu = [ -1 1 1 0];

% loop over reaction coordinate
for i = 1: size
    % calculate final mole for given reaction coordinate
    n.f = n.i  +  nu * xi(i);
    y = n.f/sum(n.f);
    Gtot(i) = sum( n.f .* ( DGoRT + log(y * P) ) );
end

plot (xi,Gtot)
xlabel('\xi_i')
ylabel('G(total)/RT')

    