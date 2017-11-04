% Parameters for the Greenland simulation

% Parameters
A0      = 1.0;
beta    = 0.005;
c       = 2e6;
CR      = 5e5;
d0      = 1545;
rc      = 8e5;
% h_Eq is set in oerlemansClassoc.m
% http://www.arctic.noaa.gov/reportcard/greenland_ice_sheet.html
hE0     = 1545;     % See oerlemansClassic.m
f       = 0.5;      % 0.5 in the fortran code - 1.0 in the article
%s       = 0.001;
s = d0 / rc;
mu0     = 8;
%mu      = mu0*c*s^2; % Equation (4) in Oerlemans (2003) REMOVE 27Sep2016
mu      = mu0 + c*s^2; % Equation (4) in Oerlemans (2003) ADD 27Sep2016
rhoi    = 900.0;
rhow    = 1025.0;
rhom    = 3500.0;
rgr     = rc;
% This gives a steady state volume of about 7 meters for T=0 in the
% non-noisy case. This corresponds to the GrIS (IPCC, WGI, Chap 13)
Tbar    = 5.8; 

% The desert effect
Afunc = @(A0, R, CR) A0.*exp(-R./CR);
%Afunc = @(A0, R, CR) A0;

% These are not needed in our simulation
hEA     = nan;
P       = nan;
xfactor = nan;

% Struct-ise the parameters we just loaded
par = struct('A0', A0, 'Afunc', Afunc, 'beta', beta, ... %'c', c,...
    'CR', CR, 'd0', d0, ...
    'f', f, 'hE0', hE0, 'hEA', hEA, 'mu', mu, 'P', P, 'rhoi', rhoi, ...
    'rhow', rhow, 'rhom', rhom, 'rc', rc, 'rgr', rgr, 's', s, ...
    'Tbar', Tbar, 'xfactor', xfactor);

% The noise
%muval       = 0.0;
%sigmaval    = 0.8576;

% Model output is in m3 ice. We want to convert to meters sea level
% http://ngdc.noaa.gov/mgg/global/etopo1_ocean_volumes.html
% Ocean covers this many square kilometers:
Ocean_km2 = 3.619e8;
% ... hence this many suare meters
Ocean_m2 = Ocean_km2 * 1e6;
% We calcuate how many meters one m3 of ice makes the ocean rise
rhoi    = 900.0;    % Density of ice
rhow    = 1025.0;   % Density of sea water
% One X ice correponds to this many X seawater
SW_From_Ice = rhoi / rhow;


% Metres sea level from one Cubic Metre Ice
OceanSurf = SW_From_Ice / Ocean_m2;