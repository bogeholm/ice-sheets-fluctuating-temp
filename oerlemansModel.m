% Integrate Equation (13) from  J. Oerlemans (2003)
% Translation of Hans Oerlemans' code from FORTRAN to MATLAB
% Troels Mikkelsen, bogeholm@nbi.ku.dk
% January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target: dRdt = f(t, R, T; parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dRdt, Btot] = oerlemansModel(time, R, T, par)


xfactor = par.xfactor;
%pi      = 3.1416
%const   = 0.75/pi;
f       = par.f;
%c       = par.c;
rhoi    = par.rhoi;
rhow    = par.rhow;
rhom    = par.rhom;
s       = par.s;
A0      = par.A0;
hE0     = par.hE0;
% Changing Tbar corresponds to changing the height of the equilibrium
% line. See first "if-else-end"-loop below.
Tbar    = par.Tbar;
d0      = par.d0;
beta    = par.beta;
rc      = par.rc;
rgr     = par.rgr;
hEA     = par.hEA;
P       = par.P;
CR      = par.CR;
Afunc   = par.Afunc;

% Calculated "once" - same from call to call
eps1    = rhoi/(rhom - rhoi);
eps2    = -rhow/(rhom - rhoi);
delta   = rhow/rhoi;
xmu     = par.mu;%8.0 + c*s^2;
sxmu    = sqrt(xmu); 


% Calculated from call to call
% Check whether we use T or P
if isnan(T)
    % This option is relevant to generate data comparable to Figure (6)
    % In Oerlemans (2003)
    hEq = d0 - hEA*sin(2.0*pi*time/P) + xfactor;
else
    % Minimal glacier models, Johannes Oerlemans, p. 8:
    % "A typical value for ?a is 0.0065 K m-1, so the change in 
    % equilibrium-line altitude per degree warming would be 154 m."
    hEq = hE0 + (T - Tbar)*1000/6.5;
end

% Accumulation rate and the desert effect
% A = A0*exp(-R/CR) or A = A0? Page 446, Equation (20) of Oerlemans (2003)
A = Afunc(A0, R, CR);

% ha = height of runoff line 
% Called h_R in the article; Equation 15
ha      = hEq + A/beta;
% he = height of the bedrock where the ice sheet ends
he      = d0 - s*R;
% Location where the runoff line intersects the ice sheet surface
% Called r_R in the article; Equation (17) 
ra      = R - (ha - he)^2/xmu;  

% Check if the ice sheet extends into the sea, ie if R > rc
% If so, use Equation (7) in the article to define the grounding line
if R > rc       % (R.gt.rc)
    rgr = R - he^2/xmu;
end

% If the radial position of the runoff line is larger than for the grounding
% line, set the runoff coordinate equal to the grounding line coordinate.
if ra > rgr     % (ra.gt.rgr)
    ra = rgr;
end

% If the height of the runoff line is less than the height of the
% the bedrock where the ice sheet ends, set the radial coordinate of the 
% runoff line (ra) to the radius of the ice sheet. This implies the same 
% accumulation rate over the entire ice sheet.
if ha < he %(ha.lt.he) 
    ra = R;
end

% Continental ice sheet?
if R <= rc
    % Avoid pathological behavior, pt. 1
    if ra < 0
        ra = 0;
    end
    % Avoid pathological behavior, pt. 2
    if R < 1    % (R.lt.1.) R=1.
        R = 1.0;
    end
    
    Btot =pi*A*R^2 ...
        -pi*beta*(ha - he)*(R^2 - ra^2) ...
        +4*pi*beta*sxmu/5.*(R - ra)^2.5 ...
        -4*pi*beta*sxmu/3.*R*(R - ra)^1.5;
end

% Ice sheet extends into the sea
if R > rc
    Btot = pi*A*rgr^2 ...
        -pi*beta*(ha - he)*(rgr^2 - ra^2) ...
        +4*pi*beta*sxmu/5.*((R - ra)^2.5-(R - rgr)^2.5) ...
        -4*pi*beta*sxmu/3.0*(R*(R - ra)^1.5 - R*(R - rgr)^1.5) ...
        -2.0*pi*rgr*delta*f*(s*rgr - d0)^2;
        % R is very close
end

% Total mass balance; Equation (13) in the article.
% This term is alway on ...
fac = pi*(1.0 + eps1)*(1.333*sxmu*R^1.5 - s*R^2);
% ... this term is only for a calving ice sheet
if R > rc
    fac = fac + 2.*eps2*(pi*s*R^2-d0*R);
end

dRdt = Btot/fac;    % Return value

end