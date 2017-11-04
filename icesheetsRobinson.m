%% Apply effect of fluctuating temperature to Robinson et al.'s results
%
% Troels B. Mikkelsen - bogeholm@nbi.ku.dk
% 2015 - 2017


% -------------------------------------------------------------------------
clearvars; close all; clc; format compact
rng('default')
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
tic
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Add subfolders in userpath
addpath(genpath(userpath));
% -------------------------------------------------------------------------


% ------------- Common setup ----------------------------------------------
run('icesheetsSetup')
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Save?
% pdf's are used in the article, png's in the README
%save_pdf = true;
%save_png = true;
save_pdf = false;
save_png = false;
% -------------------------------------------------------------------------




% ------------- Variance of AR(1) process ---------------------------------
%https://en.wikipedia.org/wiki/Autoregressive_model#...
% Example:_An_AR.281.29_process
%
% x(t+1) = a1*x(t) + s*eta => var(xt) = s^2/(1 - a1^2)
arvar_func = @(a1, s) s^2 / (1 - a1^2);
% -------------------------------------------------------------------------

% ------------- m. SLE per Gt ice (see Supplementing Information) ---------
m_SLE_pr_Gt = 2.6958e-6;
mm_SLE_pr_Gt = 2.6958e-3;
% -------------------------------------------------------------------------

% ------------- From greenlandSigma.m -------------------------------------
try
    ar1results = load([datapath, 'ar1results.mat']);
catch FE
    disp('No file ar1results.mat - run greenlandTemperature.m')
    % Can't proceed without these results
    rethrow(FE)
end
%display(ar1results.estar1)
%display(ar1results.estvar)
estar1 = ar1results.estar1;
estvar = ar1results.estvar;
arvar = arvar_func(estar1, sqrt(estvar));
sigmaval = sqrt(arvar);
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Load the parameters for this run
run('oerlemansParam')
% -------------------------------------------------------------------------



%% Load Robinsons data
datafile = [robinsonpath, 'Robinson2012_transient[December2015].nc'];
ncid = netcdf.open(datafile, 'NOWRITE');
[ndims, nvars, natts, unlimdimID] = netcdf.inq(ncid);

% Display variable names
for ii = 0:nvars-1
    [varname, vartype, dimids, natts] = netcdf.inqVar(ncid,ii);
    display(varname)
    display(natts)
end

%netcdf.close(ncid)
%netcdf.inqAttName(ncid, 2,0)
smbid           = netcdf.inqVarID(ncid, 'smb');
smbunits        = netcdf.getAtt(ncid, smbid, 'units'); 


% Query time units
timeid           = netcdf.inqVarID(ncid, 'time');
timeunits        = netcdf.getAtt(ncid, timeid, 'units'); 
% Time axis in kilo years
time_rob        = ncread(datafile, 'time');
% The results
smb_rob         = ncread(datafile, 'smb');
% -------------------------------------------------------------------------
% Robinson data in mm SLE / yr
smb_rob = smb_rob*mm_SLE_pr_Gt;
% -------------------------------------------------------------------------
Vtot_rob        = ncread(datafile, 'Vtot');

vol_id          = netcdf.inqVarID(ncid, 'Vtot');
vol_unit        = netcdf.getAtt(ncid, vol_id, 'units'); 


% Parameters
T_warming_rob   = ncread(datafile, 'T_warming');
ppfac_rob       = ncread(datafile, 'ppfac');
itm_c_rob       = ncread(datafile, 'itm_c');
dT_jja_rob      = ncread(datafile, 'dT_jja');



% There are numel(T_warming_rob) * numel(ppfac_rob) * numel(itm_c_rob) = 
% 11 x 9 x 11 = 1089 simulations.

% (1x11), NOT in order:
T_warming = [7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 3.5, 4.5, 2.5, 6.5, 5.5];
% (1x9), looks like: ppfac = [-0.07, -0.06, ..., 0.00, 0.01];
% The ppfac parameters are close but not exact to two decimal points
ppfac = unique(ppfac_rob);
% (1x11), looks like = [-60, ... -50];
itm_c = unique(itm_c_rob);

% Since the temperatures were ramped up the first 100 years, we start from
% t >= 100y; this corresponds to index 21 in the time vector.
% We give the ice sheet 100 years more to stabilize, corresponding to 
% index 41 or 200 years
first = 41;

% Easy mapping from (T, params) to simulation number
simno = @(x, y, z) simnumber(x, y, z, T_warming, ppfac, itm_c);

% The temperatures are not in order
[Tsorted, Tperm] = sort(T_warming);  
% Number of parameters
nitm = numel(itm_c);
nppf = numel(ppfac);
nT = numel(Tsorted);


%% Plot all V(t)
fig911 = figure(911); hold on; box on; figset(fig911)

%cls = lines(numel(Tsorted));
cls = redblue(numel(Tsorted));

% We want volume in meters sea level equivalent (SLE)
% SLE = OceanSurf * V, where V is in m3
%
% Vtot_rob has units 1e6 km3 = 1e15 m3
rob_sle = 1e15*OceanSurf;

for idi = 1:nitm
    for idj = 1:nppf
        for idk = 1:numel(Tsorted)
            %N = simno(Tsorted(idk), itm_c(idi), ppfac(idj));
            N = simno(Tsorted(idk), itm_c(idi), ppfac(idj));
            % NaN / collapse is set to -9999
            vt = Vtot_rob(N, :);
            vt = vt(vt ~= -9999);
            tax = time_rob(vt ~= -9999);
            plot(tax, rob_sle*vt, 'Color', cls(idk, :))
        end
    end
end

xlim([0.0 2]);
%ylim([0.0 4.0]);

xl = xlabel('Time (ka)'); textset(xl)
yl = ylabel('Volume (m. SLE)'); textset(yl)

line(0.2*[1, 1], ylim, 'Color', 'k')

%% Plot all V(t) - separate the parameters
fig117 = figure(117); hold on; box on; figset(fig117)

% disp(nitm); % 11
% disp(nppf); % 9

%cls = lines(numel(Tsorted));
cls = lines(nT);

% Consider a miximum temperature? Set to inf for not
maxtemp = 4;
numtemps = sum(Tsorted <= maxtemp);

for idx = 1:nppf
    subplot(330 + idx); hold on;
    for idy = 1:4%nitm
        for idz = 1:numtemps
            % Get the simulation number
            N = simno(Tsorted(idz), itm_c(idy), ppfac(idx));
            % NaN / collapse is set to -9999
            vt = Vtot_rob(N, :);
            vt = vt(vt ~= -9999);
            tax = time_rob(vt ~= -9999);
            plot(tax, vt, 'Color', cls(idy, :))
        end
        xlim([0.0 0.4]);
        ylim([3.3 3.8]);
        line(0.2*[1, 1], ylim, 'Color', 'k', 'Linewidth', 0.5)
    end
end


[~, xl] = suplabel('Time [kyr]', 'x'); textset(xl)
[~, yl] = suplabel('Volume, V [$10^6$ km$^3$]', 'y'); textset(yl)

%% Find out how much the V(T) varies at t_fit for each parameter pair
volatfit = cell(nppf, nitm);


% Get all volumes at fitting time
for idx = 1:nppf
    for idy = 1:nitm
        % Save results here
        vols = zeros(numel(Tsorted), 1);
        for idz = 1:numel(Tsorted)
            % Get the simulation number
            N = simno(Tsorted(idz), itm_c(idy), ppfac(idx));
            % NaN / collapse is set to -9999
            vt = Vtot_rob(N, first);
            vols(idz, 1) = vt;
        end
        volatfit{idx, idy} = vols;
    end
end

% Get the max difference at fitting time
fitdiffs = nan(nppf, nitm);
fitmeans = nan(nppf, nitm);

% Get all volumes at fitting time
for idx = 1:nppf
    for idy = 1:nitm
        % Save results here
        vols = volatfit{idx, idy};
        % Save the difference, means and std's
        fitdiffs(idx, idy) = max(vols) - min(vols);
        fitmeans(idx, idy) = mean(vols);
    end
end


% Histograms of the above
edges = 4:0.5:10;
fig118 = figure(118); figset(fig118)
hhand = histogram(fitdiffs(:) ./ fitmeans(:)*100, edges);

disp('*** Histogram - all values - number of values')
disp(numel(fitdiffs))

xl = xlabel('(Max. diff. in V) / (mean V) in percent'); textset(xl)
yl = ylabel('Number of parameter pairs'); textset(yl)



%% Find out how much the V(T) varies at t_fit for each parameter pair
%  Only consider warmings of maxtemp or less
volatfit_var = cell(nppf, nitm);


% Get all volumes at fitting time
for idx = 1:nppf
    for idy = 1:nitm
        % Save results here
        vols = zeros(numtemps, 1);
        for idz = 1:numtemps
            % Get the simulation number
            N = simno(Tsorted(idz), itm_c(idy), ppfac(idx));
            % NaN / collapse is set to -9999
            vt = Vtot_rob(N, first);
            vols(idz, 1) = vt;
        end
        volatfit_var{idx, idy} = vols;
    end
end

% Get the max difference at fitting time
fitdiffs = nan(nppf, nitm);
fitmeans = nan(nppf, nitm);

% Get all volumes at fitting time
for idx = 1:nppf
    for idy = 1:nitm
        % Save results here
        vols = volatfit_var{idx, idy};
        % Save the difference, means and std's
        fitdiffs(idx, idy) = max(vols) - min(vols);
        fitmeans(idx, idy) = mean(vols);
    end
end


% Histograms of the above
edges = 0:0.25:3;
fig119 = figure(119); figset(fig119)
hhand = histogram(fitdiffs(:) ./ fitmeans(:)*100, edges);

disp('*** Histogram - some values - number of values')
disp(numel(fitdiffs))

xl = xlabel('(Max. diff. in V) / (mean V) in percent'); textset(xl)
yl = ylabel('Number of parameter pairs'); textset(yl)



%% Plot all T(t) for dt = 5 years and t > t(first)
fig112 = figure(112); hold on; box on; figset(fig112)

cls = redblue(numel(Tsorted));

for idi = 1:nitm
    for idj = 1:nppf
        for idk = 1:numel(Tsorted)
            N = simno(Tsorted(idk), itm_c(idi), ppfac(idj));
            % NaN / collapse is set to -9999
            t = dT_jja_rob(N, :);
            t = t(t ~= -9999);
            tax = time_rob(t ~= -9999);
            plot(tax, t, 'Color', cls(idk, :))
        end
    end
end

xlim([0 1]);
ylim([0 8]);

xl = xlabel('Time (ka)'); textset(xl)
yl = ylabel('JJA Mean Temp. Anomaly ($^{\circ}$C)'); textset(yl)

%% Fit 3rd degree polynomial model to all SMBs
fprintf('Fitting polynomials...\n');

% Polynomial model, 3rd degree
ftypePoly3 = fittype('A*T^3 + B*T^2 + C*T + D', ...
    'independent', {'T'}, ...
    'coefficients', {'A', 'B', 'C', 'D'});
optPoly3 = fitoptions(ftypePoly3);
optPoly3.Startpoint = [1, 1, 1, 1];

tthis = first;

% Store SMB's here
SMBmat= zeros(nitm, nppf, nT);

% Store A, B and C here
Amat = zeros(nitm, nppf);
Bmat = zeros(nitm, nppf);
Cmat = zeros(nitm, nppf);
Dmat = zeros(nitm, nppf);

for ii = 1:nitm
    for jj = 1:nppf
        SMBfT = zeros(size(Tsorted));
        for kk = 1:numel(Tsorted)
            N = simno(Tsorted(kk), itm_c(ii), ppfac(jj));
            SMBfT(kk) = smb_rob(N, tthis);
        end
        SMBmat(ii, jj, :) = SMBfT;
        smbPoly = fit(Tsorted', SMBfT', ftypePoly3, optPoly3);
        % Hide the coeffs
        Amat(ii, jj) = smbPoly.A;
        Bmat(ii, jj) = smbPoly.B;
        Cmat(ii, jj) = smbPoly.C;
        Dmat(ii, jj) = smbPoly.D;
    end
end

avec = reshape(Amat, 1, nitm*nppf);
bvec = reshape(Bmat, 1, nitm*nppf);
cvec = reshape(Cmat, 1, nitm*nppf);
dvec = reshape(Dmat, 1, nitm*nppf);



%{
%% Check to see we did it right
Ntry = 2;
ivec = datasample(1:nitm, Ntry, 'Replace', false);
jvec = datasample(1:nppf, Ntry, 'Replace', false);
colors = lines(Ntry^2);

fig001 = figure(001); hold on; box on; figset(fig001)
for idi = 1:Ntry
    ii = ivec(idi);
    for idj = 1:Ntry
        jj = jvec(idj);
        hold on
        col = colors(idi + (idj-1)*Ntry, :);
        scatter(Tsorted, reshape(SMBmat(ii, jj, :), 1, nT), [], col,...
            'filled')
        
        func =  @(T) Amat(ii, jj)*T.^3 + Bmat(ii, jj)*T.^2 + ...
            Cmat(ii, jj)*T + Dmat(ii, jj);
        plot(Tsorted, func(Tsorted), 'Color', col)
        
        %{
        funcTT =  @(T) 6*Amat(ii, jj)*T + 2*Bmat(ii, jj);
        
        plot(Tsorted, func(Tsorted) + funcTT(Tsorted),...
            'Color', col, 'Linestyle', ':')
        %}
    end
end


% Labels, fig color
xl = xlabel('Warming T [$^{\circ}$C]');
yl = ylabel('SMB [mm SLE yr$^{-1}$]');
%tt = title('Polynomial fit to SMB(T)'); textset(tt);
textset(xl); textset(yl);

%close(678) % Only a check!


%}

%% Showcase the minimization objective
%  Choose an arbitrary parameter set from avec, bvec, ...
rng('default');
idx = datasample(1:nitm*nppf, 1, 'Replace', false);
%idx = 19;
% Choose a temperature - this one sits nicely in the middle of the figure
tplot = 4;
% T is for plotting - increasuing resolution will increase correctness
% of minimization
T = linspace(min(Tsorted), max(Tsorted), 1000);

% The mass balance function "f of T"
smb =  @(T) avec(idx)*T.^3 + bvec(idx)*T.^2 + ...
    cvec(idx)*T + dvec(idx);
% The second derivative wrt temperature "f_TT of T"
smbacc =  @(T) 6*avec(idx)*T + 2*bvec(idx);
% Add the two functions, "function (two) of T "
smbtot = @(T) smb(T) + 0.5*sqrt(arvar)*smbacc(T);

% Plot limits should be adjusted accordingly (only visual effetct).
manual_limits = false;
if manual_limits
    xlimits = [3 5];
    ylimits = [-2.5 -0.2];
else
    xlimits = [3.35 4.5];
    ylimits = [smb(tplot)-0.45, smb(tplot)+0.25];
end

% Get x, y distances
xdist = xlimits(2) - xlimits(1);
ydist = ylimits(2) - ylimits(1);

% Plot f, f2
fig002 = figure(002); hold on; box on; figset(fig002)
pf0 = plot(T, smb(T), 'Color', blue, 'linewidth', 1.8);
pf2 = plot(T, smbtot(T), 'Color', blue, 'linewidth', 1.8, ...
    'Linestyle', ':');
xlim(xlimits);
ylim(ylimits);



% -------------------------------------------------------------------------
% Simple minimization procedure
[~, indmin] = min(abs(smbtot(T) - smb(tplot)));
t_other = T(indmin);
% Accuracy increases with resolution of T
% -------------------------------------------------------------------------


% Print some results
del_t = tplot - t_other;
del_smb = smb(tplot) - smbtot(tplot);

fprintf('g(T_0) = %.2f\n',    smb(tplot));
fprintf('g*(T_0) = %.2f\n',   smbtot(tplot));
fprintf('T* = %.2f\n',        t_other);
fprintf('g*(T*) = %.2f\n',    smbtot(t_other));
fprintf('delta T = %.2f\n',   del_t);
fprintf('delta smb = %.2f\n', del_smb);



% Scatter the results of minimization
markersize = 100;
% (T_0, SMB(T_0))
scf0 =  scatter(tplot, smb(tplot), markersize, 'o', ...
    'markerfacecolor', green, 'markeredgecolor', green);
% (T_0, SMB^*(T_0))
scf2 = scatter(tplot, smbtot(tplot), markersize, 's', ...
    'markerfacecolor', purple, 'markeredgecolor', purple);
% (T^*, SMB^*(T^*))
scfn = scatter(t_other, smbtot(t_other), markersize, 'h', ...
    'markerfacecolor', red, 'markeredgecolor', red);
  



% Showcase these concepts on the figure
% Vertical arrow
xv = tplot*[1 1];
yv = [smbtot(tplot) smb(tplot) ];
% Normalized figure coordinates
[xnv, ynv] = ds2nfu(xv, yv);
haxv = annotation('arrow', xnv, ynv);
% Horizontal arrow
xh = [t_other tplot];
yh = [smbtot(t_other), smb(tplot)];
% Normalized figure coordinates
[xnh, ynh] = ds2nfu(xh, yh);
haxh = annotation('arrow', xnh, ynh);


% (x, y)-limits and distance
xl = xlim; 
xd = xl(2) - xl(1);
yl = ylim; 
yd = yl(2) - yl(1);


% Textarrow specifying delta SMB
xa1 = [tplot+0.2*xd tplot ];
% Midpoint of vertical arrow
yp = 0.5*(smb(tplot) + smbtot(tplot));
% y-coordinate
ya1 = [yp - 0.0*yd yp];
% Normalized figure coordinates
[nxa1, nya1] = ds2nfu(xa1, ya1);
% The annotation itself
sa1 = '$\Delta$SMB$ > 0$';
htex1 = annotation('textarrow', nxa1, nya1, 'String', sa1); textset(htex1)


% Textarrow specifying delta T
% Midpoint of vertical arrow
xp = 0.5*(tplot + t_other);
xa2 = [xp + 0.13*xd xp];
% y-coordinate
ya2 = [0.2*yd + smbtot(t_other) smbtot(t_other)];
% Normalized figure coordinates
[nxa2, nya2] = ds2nfu(xa2, ya2);
% The annotation itself
sa2 = '$\Delta T > 0$';
htex2 = annotation('textarrow', nxa2, nya2, 'String', sa2); textset(htex2)

% -------------------------------------------------------------------------
% Legend and text, CASE 1: f(T_0) + f_TT(T_0) ...
% -------------------------------------------------------------------------
%{
% Legend strings
legstrs = {'$\tilde{f}(T)$ := SMB$(T)$', ...
    '$\tilde{f}(T) + \sigma_T^2/2 \cdot \tilde{f}_{TT}(T)$', ...
    '$\tilde{f}(T_0)$', ...
    '$\tilde{f}(T_0) + \sigma_T^2/2 \cdot \tilde{f}_{TT}(T_0)$', ...
    '$\tilde{f}(T_0 + \Delta T) + \sigma_T^2/2 \cdot \tilde{f}_{TT}(T_0 + \Delta T)$'};
    %'$\tilde{f}(\hat{T}) + \sigma_T^2/2 \cdot \tilde{f}_{TT}(\hat{T})$'};
    
% Legend
l1 = legend([pf0 pf2 scf0 scf2 scfn], legstrs, ...
    'location', 'southwest'); legset(l1)

% Put text at the various points
tp0 = text(tplot+0.05*xdist, smb(tplot), ...
    'SMB$(T_0)$', 'color', green); textset(tp0)

tp1 = text(tplot-0.3*xdist, smbtot(tplot)-0.0*ydist, ...
    'SMB$(T_0) + \Delta$SMB$(T_0)$', ...
    'color', purple); textset(tp1)

% LaTeX newline in MATLAB is not trivial
tp2_0 = text(t_other-0.25*xdist, smbtot(t_other)-0.0*ydist,...
      'SMB$(T_0 + \Delta T)$', ...
      'color', red); textset(tp2_0);

% LaTeX newline in MATLAB is not trivial
tp2_1 = text(t_other-0.25*xdist, smbtot(t_other)-0.05*ydist,...
      '$+\Delta $SMB$ (T_0 + \Delta T)$', ...
      'color', red); textset(tp2_1);
%}
% -------------------------------------------------------------------------
% End case 1
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Legend and text, CASE 2: [f + f_TT](T_0) ...
% Legend strings
legstrs = {'$\tilde{f}(T)$ := SMB$(T)$', ...
    '$\tilde{f}(T_0)$', ...
    '$\left[\tilde{f} + \sigma_T^2/2 \cdot \tilde{f}_{TT} \right](T)$', ...
    '$\left[\tilde{f} + \sigma_T^2/2 \cdot \tilde{f}_{TT}\right](T_0)$', ...
    '$\left[\tilde{f} + \sigma_T^2/2 \cdot \tilde{f}_{TT}\right](T_0 - \Delta T)$'};
    %'$\tilde{f}(\hat{T}) + \sigma_T^2/2 \times \tilde{f}_{TT}(\hat{T})$'};
    
% Legend
l1 = legend([pf0 scf0 pf2 scf2 scfn], legstrs, ...
    'location', 'southwest'); legset(l1)

% Put text at the various points
tp0 = text(tplot+0.05*xdist, smb(tplot), ...
    'SMB$(T_0)$', 'color', green); textset(tp0)

tp1 = text(tplot-0.3*xdist, smbtot(tplot)-0.0*ydist, ...
    '$[$SMB$ - \Delta$SMB$](T_0)$', ...
    'color', purple); textset(tp1)

% LaTeX newline in MATLAB is not trivial
tp2_0 = text(t_other-0.4*xdist, smbtot(t_other)-0.0*ydist,...
      '$[$SMB$ - \Delta$SMB$](T_0 - \Delta T)$', ...
      'color', red); textset(tp2_0);
% -------------------------------------------------------------------------
% End case 2
% -------------------------------------------------------------------------


% Pretty
% Larger font since this will be a sub panel
xl = xlabel('Warming, T [$^{\circ}$C]'); %textset(xl)
set(xl, 'Fontsize', 22, 'Interpreter', 'Latex');

yl = ylabel('$dV/dt$ [mm SLE yr$^{-1}$]'); %textset(yl)
set(yl, 'Fontsize', 22, 'Interpreter', 'Latex');

%print(fig002, '/Users/bogeholm/Desktop/tezt.png', '-dpng', '-r400')



%% Showcase the minimization procedure
%  Choose an arbitrary parameter set from avec, bvec, ...
idx = datasample(1:nitm*nppf, 1, 'Replace', false);
% T is for plotting - increasuing resolution will increase correctness
% of minimization
T = linspace(min(Tsorted), max(Tsorted), 1000);

% Extend T for the minimization procedure
dT = T(2) - T(1);
Tp = 0.5:dT:(2-dT);
T_ext = [Tp T];


% The mass balance function "f of T"
smb =  @(T) avec(idx)*T.^3 + bvec(idx)*T.^2 + ...
    cvec(idx)*T + dvec(idx);
% The second derivative wrt temperature "f_TT of T"
smbacc =  @(T) 6*avec(idx)*T + 2*bvec(idx);
% Add the two functions, "function (two) of T "
smbtot = @(T) smb(T) + 0.5*sqrt(arvar)*smbacc(T);

% Now compute the smb and smb acceleration
smbv = smb(T);
%smbv_e = smb(Text);
smbtv = smbtot(T);
smbtv_e = smbtot(T_ext);


del_t = zeros(size(T));
del_smb = zeros(size(T));

% Loop over T, do minimization
for idx = 1:numel(T)
    t_here = T(idx);
    smb_here = smbv(idx);
    
    % Simple minimization procedure
    [~, indmin] = min(abs(smbtv_e - smb_here));
    t_other = T_ext(indmin);
    % Accuracy increases with resolution of T
    del_t_val   = t_here - t_other;
    del_smb_val = smb_here - smbtv(idx);
    
    del_t(idx) = del_t_val;
    del_smb(idx) = del_smb_val;
    

end


%{
% This is just a figure of the SMB here
fig003 = figure(003); hold on; box on; figset(fig003)
pf0 = plot(T, smb(T), 'Color', blue, 'linewidth', 1.8);
pf2 = plot(T, smbtot(T), 'Color', blue, 'linewidth', 1.8, 'Linestyle', ':');
% Pretty
xl = xlabel('Warming, T [$^{\circ}$C]'); textset(xl)
yl = ylabel('SMB [mm SLE yr$^{-1}$]'); textset(yl)
%}

% Plot the delta T and delta SMB
fig004 = figure(004); hold on; box on; figset(fig004)
% First, Delta T
subplot(1, 2, 1)
plot(T, del_t)
yl = ylabel('$\Delta T$'); textset(yl)
% Now, delta SMB
subplot(1, 2, 2)
plot(T, del_smb)
yl = ylabel('$\Delta$ SMB [mm SLE yr$^{-1}$]'); textset(yl)
% Common xlabel
[~, xl] = suplabel('Warming, T [$^{\circ}$C]', 'x'); textset(xl)





%% Run minimization procedure for ALL parameters
fprintf('Minimizing to find T*...\n');

% Adding one ensures that {2,3,4,...} is divisible by dT 
T = linspace(min(Tsorted), max(Tsorted), 2001); 
% Extend T for the minimization procedure
dT = T(2) - T(1);
Tp = 0.5:dT:(2-dT);
T_ext = [Tp T];

% Store results
Del_T = zeros(numel(T), nitm*nppf);
Del_SMB = zeros(numel(T), nitm*nppf);

for idp = 1:nitm*nppf
    % The mass balance function "f of T"
    smb =  @(T) avec(idp)*T.^3 + bvec(idp)*T.^2 + ...
        cvec(idp)*T + dvec(idp);
    % The second derivative wrt temperature "f_TT of T"
    smbacc =  @(T) 6*avec(idp)*T + 2*bvec(idp);
    % Add the two functions, "function (two) of T "
    smbtot = @(T) smb(T) + 0.5*sqrt(arvar)*smbacc(T);
    
    
    % Now compute the smb and smb acceleration
    smbv = smb(T);
    %smbv_e = smb(Text);
    smbtv = smbtot(T);
    smbtv_e = smbtot(T_ext);
    
    
    for idt = 1:numel(T)
        t_here = T(idt);
        smb_here = smbv(idt);

        % Simple minimization procedure
        [~, indmin] = min(abs(smbtv_e - smb_here));
        t_other = T_ext(indmin);
        % Accuracy increases with resolution of T
        del_t_val   = t_here - t_other;
        del_smb_val = smb_here - smbtv(idt);

        Del_T(idt, idp) = del_t_val;
        Del_SMB(idt, idp) = del_smb_val;

    end
end

% These figures will be plotted better below
%{
% Computed Delta T
fig005 = figure(005); hold on; box on; figset(fig005)
plot(T, Del_T)
% Pretty
xl = xlabel('Warming, T [$^{\circ}$C]'); textset(xl)
yl = ylabel('$\Delta$T [$^{\circ}$C]'); textset(yl)


% Computed Delta SMB
fig006 = figure(006); hold on; box on; figset(fig006)
plot(T, Del_SMB)
% Pretty
xl = xlabel('Warming, T [$^{\circ}$C]'); textset(xl)
yl = ylabel('$\Delta$SMB [mm SLE yr$^{-1}$]'); textset(yl)
%}

%% Compute credible intervals of Delta T and Delta SMB
fprintf('Computing credible intervals...\n');

% Credibility level
alpha = 0.05;
% Number of points for ksdensity()
kspoints = 500;

% Store results for Delta T here
Del_T_likely = zeros(numel(T), 1);
Del_T_min = zeros(numel(T), 1);
Del_T_max = zeros(numel(T), 1);

% Loop over temperatures - do ksdensity for each T
for idt = 1:numel(T)
    
    % All values of delta T for this temperature
    del_t_here = Del_T(idt, :);
    
    [fks_t, xks_t]      = ksdensity(del_t_here, 'npoints', kspoints);
    
    % Maximum
    [~, ilikely_t]      = max(fks_t);
    
    % Construct credible intervals
    ct_t = cumtrapz(xks_t, fks_t);
    imin_t = find(ct_t < alpha/2, 1, 'last');
    imax_t = find(ct_t > 1 - alpha/2, 1, 'first') + 1;
    
    % Get values using indices
    likely_del_t = xks_t(ilikely_t);
    min_del_t = xks_t(imin_t);
    max_del_t = xks_t(imax_t);
    
    Del_T_likely(idt) = likely_del_t;
    Del_T_min(idt) = min_del_t;
    Del_T_max(idt) = max_del_t;
    
end

% Same as above, for Delta SMB
Del_SMB_likely = zeros(numel(T), 1);
Del_SMB_min = zeros(numel(T), 1);
Del_SMB_max = zeros(numel(T), 1);

for idt = 1:numel(T)
    % All values of delta T for this temperature
    del_smb_here = Del_SMB(idt, :);
    
    [fks_smb, xks_smb]      = ksdensity(del_smb_here, 'npoints', kspoints);
    
    % Maximum
    [~, ilikely_smb]      = max(fks_smb);
    
    % Construct credible intervals
    ct_smb = cumtrapz(xks_smb, fks_smb);
    imin_smb = find(ct_smb < alpha/2, 1, 'last');
    imax_smb = find(ct_smb > 1 - alpha/2, 1, 'first') + 1;
    
    % Get values using indices
    likely_del_smb = xks_smb(ilikely_smb);
    min_del_smb = xks_smb(imin_smb);
    max_del_smb = xks_smb(imax_smb);
    
    Del_SMB_likely(idt) = likely_del_smb;
    Del_SMB_min(idt) = min_del_smb;
    Del_SMB_max(idt) = max_del_smb;
    
end


%% Output table of delta SMB, delta T
% Do it for these temperatures
Tinterest = [2, 3, 4, 5, 6, 7];
indices = zeros(size(Tinterest));


for idx = 1:numel(indices)
    indhere = find(T == Tinterest(idx));
    assert(~isempty(indhere), 'Could not find that - try to change dT')
    indices(idx) = indhere;
end


Del_T_table = zeros(numel(indices), 4);
Del_SMB_table = zeros(numel(indices), 4);

for idx = 1:numel(indices)
    index = indices(idx);
    % Make table for delta T
    Del_T_table(idx, 1) = T(index);
    Del_T_table(idx, 2) = Del_T_likely(index);
    Del_T_table(idx, 3) = Del_T_min(index);
    Del_T_table(idx, 4) = Del_T_max(index);

    % Make table for delta T
    Del_SMB_table(idx, 1) = T(index);
    Del_SMB_table(idx, 2) = Del_SMB_likely(index);
    Del_SMB_table(idx, 3) = Del_SMB_min(index);
    Del_SMB_table(idx, 4) = Del_SMB_max(index);
    
    
end


fprintf(' -------------- Delta T -------------- \n')
fprintf('\t   T \t   Delta T \t   Lower \t    Upper \n')
format bank % two decimal points; 'currency' format
display(Del_T_table);
fprintf('\n')
fprintf(' -------------- Delta SMB ------------ \n')
fprintf('\t   T \t   Delta SMB \t   Lower \t    Upper \n')
display(Del_SMB_table);
format short % default

%% Plot Delta SMB and Delta T with credible intervals
% Figure: Delta T
fig007 = figure(007); hold on; box on; figset(fig007)
% Fill can behave highly unintuitively (to me ;-) - check with the yellow 
% lines to be sure
XT = [T, fliplr(T)];
YT = [Del_T_min', fliplr(Del_T_max')];
%plot(T, Del_T_min, 'color', yellow, 'linewidth', 5)
%plot(T, Del_T_max, 'color', yellow, 'linewidth', 5)

facealpha = 0.7;            % Alpha for the blue region
lightgrey = 0.25*[1 1 1];   % Color for the grey lines
greywidth = 0.75;            % Linewidth for the grey lines
redwidth = 2.5;               % Linewidth for the red lines


% Fill 95% credible region
p1 = plot(T, Del_T, 'color', lightgrey, 'linewidth', greywidth);
fi1 = fill(XT, YT, 'b');
set(fi1, 'FaceColor', blue, 'FaceAlpha', facealpha, 'EdgeColor', 'none')
p2 = plot(T, Del_T_likely, 'color', violetred, 'linewidth', redwidth);

% Legend
legstrs = {'Most likely $\Delta$T',...
    'All $\Delta$T', ...
    '$95\%$ Credibility region'};
l1 = legend([p2 p1(1) fi1], legstrs); legset(l1)


% Pretty
xl = xlabel('Warming, T [$^{\circ}$C]'); textset(xl)
yl = ylabel('$\Delta$T [$^{\circ}$C]'); textset(yl)


% Figure: Delta SMB
fig008 = figure(008); hold on; box on; figset(fig008)
% Fill can behave highly unintuitively (to me ;-) - check with the yellow 
% lines to be sure. Get the apostrope INSIDE fliplr
XSMB = [T, fliplr(T)];
YSMB = [Del_SMB_min', fliplr(Del_SMB_max')];
%plot(T, Del_T_min, 'color', yellow, 'linewidth', 5)
%plot(T, Del_T_max, 'color', yellow, 'linewidth', 5)

% Fill 95% credible region
p1 = plot(T, Del_SMB, 'color', lightgrey, 'LineWidth', greywidth);
fi1 = fill(XSMB, YSMB, 'b');
set(fi1, 'FaceColor', blue, 'FaceAlpha', facealpha, 'EdgeColor', 'none')
p2 = plot(T, Del_SMB_likely, 'color', violetred, 'linewidth', redwidth);

% Legend
legstrs = {'Most likely $\Delta$SMB',...
    'All $\Delta$SMB', ...
    '$95\%$ Credibility region'};
l1 = legend([p2 p1(1) fi1], legstrs); legset(l1)

% Pretty
xl = xlabel('Warming, T [$^{\circ}$C]'); textset(xl)
yl = ylabel('$\Delta$SMB [mm SLE yr$^{-1}$]'); textset(yl)


% Figure: Delta T and Delta SMB on same plot
fig009 = figure(009); hold on; box on; figset(fig009)
subplotfont = 10;
% Delta T on top
subplot(2, 1, 1); hold on; box on
plot(T, Del_T, 'color', lightgrey, 'LineWidth', greywidth);
% Fill 95% credible region
fi1 = fill(XT, YT, 'b');
set(fi1, 'FaceColor', blue, 'EdgeColor', 'none', 'FaceAlpha', facealpha)
% Plots and lines
plot(T, Del_T_likely, 'color', violetred, 'linewidth', redwidth);
% Pretty
xlim([2 7])
yl = ylabel('$\Delta$T [$^{\circ}$C]'); %textset(yl)
set(yl, 'Fontsize', subplotfont, 'Interpreter', 'Latex');

% tick size
set(gca, 'FontSize', subplotfont);

% Delta SMB on the bottom
subplot(2, 1, 2); hold on; box on
plot(T, Del_SMB, 'color', lightgrey, 'LineWidth', greywidth);
% Fill 95% credible region
fi1 = fill(XSMB, YSMB, 'b');
set(fi1, 'FaceColor', blue, 'EdgeColor', 'none', 'FaceAlpha', facealpha)
% Plots and lines
plot(T, Del_SMB_likely, 'color', violetred, 'linewidth', redwidth);
% Pretty
xlim([2 7])
xl = xlabel('Warming, T [$^{\circ}$C]'); %textset(xl)
yl = ylabel('$\Delta$SMB [mm SLE yr$^{-1}$]'); %textset(yl)
set(yl, 'Fontsize', subplotfont, 'Interpreter', 'Latex');
set(xl, 'Fontsize', subplotfont, 'Interpreter', 'Latex');

% tick size
set(gca, 'FontSize', subplotfont);

%% Same as above - transparent!
% Figure: Delta T and Delta SMB on same plot
fig010 = figure(010); hold on; box on; figset(fig010)

% Delta T on top
subplot(2, 1, 1); hold on; box on
% Fill doesn't work with transparent background
% Fill 95% credible region
%fi1 = fill(XT, YT, 'b');
%set(fi1, 'FaceColor', blue, 'EdgeColor', 'none')
% Plots and lines
plot(T, Del_T, 'color', lightgrey, 'LineWidth', greywidth);
plot(T, Del_T_likely, 'color', violetred, 'linewidth', redwidth);
% Pretty
xlim([2 7])
yl = ylabel('$\Delta$T [$^{\circ}$C]'); textset(yl)



% Delta SMB on the bottom
subplot(2, 1, 2); hold on; box on
plot(T, Del_SMB, 'color', lightgrey, 'LineWidth', greywidth);
% Fill doesn't work with transparent background
% Fill 95% credible region
%fi1 = fill(XSMB, YSMB, 'b');
%set(fi1, 'FaceColor', blue, 'EdgeColor', 'none')
% Plots and lines
plot(T, Del_SMB_likely, 'color', violetred, 'linewidth', redwidth);
% Pretty
xlim([2 7])
xl = xlabel('Warming, T [$^{\circ}$C]'); textset(xl)
yl = ylabel('$\Delta$SMB [mm SLE yr$^{-1}$]'); textset(yl)


% pdf's
if save_pdf == true
    disp('Saving pdfs...')
    export_fig(fig002, [pdfpath, '2016gl070016-p03.pdf'])
    % "export_fig currently supports transparent patches/areas only 
    % in PNG output.
    print(fig009, [pdfpath, '2016gl070016-p04.pdf'], ...
        '-dpdf', '-r400')
    export_fig(fig112, [pdfpath, 'RobinsonTemperature.pdf'])
    export_fig(fig118, [pdfpath, 'VolumeHistogram.pdf'])
    export_fig(fig119, [pdfpath, 'VolumeHistogramMaxtemp.pdf'])
    export_fig(fig911, [pdfpath, 'RobinsonVolumeTime.pdf'])
end

% png's for the README
if save_png == true
    disp('Saving pngs...')
    print(fig002, [pngpath, '2016gl070016-p03.png'], '-dpng', '-r400')
    print(fig002, [pngpath, '2016gl070016-p03.eps'], '-depsc', '-r400')
    print(fig002, [pngpath, 'figure-two-right-panel.jpg'], '-djpeg', '-r600')
    %export_fig(fig002, [pngpath, '2016gl070016-p03.png'])
    export_fig(fig009, [pngpath, '2016gl070016-p04.png'])
    export_fig(fig118, [pngpath, 'VolumeHistogram.png'])
    export_fig(fig119, [pngpath, 'VolumeHistogramMaxtemp.png'])
end


% -------------------------------------------------------------------------
fprintf('Done.\n');
toc
% -------------------------------------------------------------------------


%% Print out data
fprintf('\n\n ******************* Computed Values *******************\n\n')

fprintf('\n-------- SMB Calculation Year ------------------------------\n')
% time_rob is in kilo years
fprintf('Year: %i\n', time_rob(first)*1000)

fprintf('\n--------------- Delta T ------------------------------------\n')
fprintf('\t   T \t   Delta T \t   Lower \t    Upper \n')
format bank % two decimal points; 'currency' format
disp(Del_T_table);

fprintf('\n--------------- Delta SMB [mm SLE/yr] ----------------------\n')
fprintf('\t   T \t   Delta SMB \t   Lower \t    Upper \n')
disp(Del_SMB_table);
format short % default

Del_SMB_Gt = Del_SMB_table;
Del_SMB_Gt(:, 2:4) = Del_SMB_Gt(:, 2:4)/mm_SLE_pr_Gt;

% Enable/disable banking format
format bank
fprintf('\n--------------- Delta SMB [Gt/yr] --------------------------\n')
fprintf('\t   T \t   Delta SMB \t   Lower \t    Upper \n')
disp(Del_SMB_Gt);
format short % default


% Delta SMB in percent of average 2003 - 2011
%  Average SMB 2003 - 2011 = -234 Gt/y (Barletta et al., 2013)
BarlettaSMB = 234; % Gt/y

format bank
fprintf('\n--------------- Delta SMB in percent of 2003 - 2011 average \n')
fprintf('\t   MLE \t   Lower \t      Upper \n')
disp(Del_SMB_Gt(2, 2:4) / BarlettaSMB * 100);
format short % default






