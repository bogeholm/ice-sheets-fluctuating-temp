function [Vcont, Vsea, Vtot] = iceVolume(R, par)

rhoi    = par.rhoi;
rhow    = par.rhow;
rhom    = par.rhom;
s       = par.s;
d0      = par.d0;
rc      = par.rc;

% Calculated "once" - same from call to call
eps1    = rhoi/(rhom - rhoi);
eps2    = -rhow/(rhom - rhoi);
mu     = par.mu;

% The continental part - Equation 9 in Oerlemans
Vcont = 8*pi*sqrt(mu)/15 * R.^(5/2) - 1/3*pi*s*R.^3;

% The sea part
Vsea = pi*(2/3*s*(R.^3 - rc^3) - d0*(R.^2 - rc^2));
% The sea part is only relevant for R > rc, otherwise formula gives
% non-sensical results
Vsea(R <= rc) = 0;

% Total volume; Oerlemans Equation 11a
Vtot = (1 + eps1)*Vcont + eps2*Vsea;