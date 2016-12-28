%% Integrate Hans Oerlemans' (2003) models with noise
%
% Troels B. Mikkelsen - bogeholm@nbi.ku.dk
% 2015 - 2016
%
%



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


% ------------- Load parameters -------------------------------------------
run('oerlemansParam')
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Save figures?
% pdf's are used in the article, png's in the README
save_pdf = true;
save_png = true;
%save_pdf = false;
%save_png = false;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%savedata = false;
savedata = true;
%thisrun = 'test';
thisrun = 'production';
% -------------------------------------------------------------------------




% ------------- Save data? ------------------------------------------------
if savedata == true
    display('*--------------> Saving data <-----------------*')
end
% -------------------------------------------------------------------------




% ------------- Parameters for integration --------------------------------
dt           = 1.0;         % 1.0 in article
timemax      = 1e5;         % 5e4 in article
% -------------------------------------------------------------------------


% ------------- Other parameters ------------------------------------------
if isequal(thisrun, 'test')
    n_sims      = 2; % Number of simulations per temperature
    %temps_show = linspace(-2, 1.5, 3)';
    % Temps for steady numerical simulations plot
    temps = [-1.5; 0; 1.5; 3];
    temps_show = temps;
    npoints      = 300;         % 
elseif isequal(thisrun, 'production')
    n_sims      = 5; % Number of simulations per temperature
    %temps_show = (-2:0.125:4)';
    % Temps for steady numerical simulations plot
    %temps = [-1.5; -0.5; 0.5; 1.5; 2.5];
    temps = [-1.5; 0; 1.5; 3];
    npoints      = 1200;         % 800 in article
    temps_show = temps;
end


n_temps     = numel(temps); % This many temperatures for OU tryouts
%n_temps_show = numel(temps_show);
% Steady state is assumed after this simulation year
steady      = 5e4;
% Average the last 'last' values
last        = floor((timemax - steady)/dt);
% -------------------------------------------------------------------------



% ------------- Variance of AR(1) process ---------------------------------
%https://en.wikipedia.org/wiki/Autoregressive_model#Example:_An_AR.281.29_process
%
% x(t+1) = a1*x(t) + s*eta => var(xt) = s^2/(1 - a1^2)
arvar_func = @(a1, s) s^2 / (1 - a1^2);
% -------------------------------------------------------------------------


% ------------- From greenlandSigma.m -------------------------------------
try
    ar1results = load([datapath, 'ar1results.mat']);
catch FE
    display('No file ar1results.mat - run greenlandTemperature2016.m')
    % Can't proceed without these results
    rethrow(FE)
end
%display(ar1results.estar1)
%display(ar1results.estvar)
estar1 = ar1results.estar1;
estvar = ar1results.estvar;
arvar = arvar_func(estar1, sqrt(estvar));
% -------------------------------------------------------------------------





%% Integrate O/03 for different steady temperatures
display('Integrating steady temperatures')
% Volume as function of radius; 
Vr = @(R) iceTotal(R, par);

Rstart = 5e5;   % Initial volume; inspect results to ensure sanity
timeax = (0:dt:timemax)';
Rresults = nan(numel(timeax), n_temps);
steps = numel(timeax);

% Run the forward Euler on O/03
for idT = 1:n_temps
    Tsim = temps(idT);
    Rval    = Rstart; % initial conditions
    time    = -dt;
    ii      = 0;
      % Do forward euler
     while time < timemax
        time    = time + dt;
        ii      = ii + 1;
        [dRdt, dVdt] = oerlemansModel(time, Rval, Tsim, par);
        Rval = Rval + dt*dRdt;
        % Store
        Rresults(ii, idT) = Rval;
     end % end while
end

% O/03 works in radius - we convert
Vresults = OceanSurf*Vr(Rresults);





%% Compute mean volume of steady temperature simulations
meanV = zeros(n_temps, 1);
for idx = 1:n_temps
    meanV(idx) = mean(Vresults((steps-last):end, idx));
end





%% Integrate O/03 for different AR(1) fluctuating temperatures 
display('Integrating OU fluctuating temperatures')

% Run the forward Euler on O/03 - with OU noise
VresAR = cell(n_temps, n_sims);
TempsAR = cell(n_temps, n_sims);

calstds = zeros(n_temps, n_sims);

for idT = 1:n_temps % First, iterate over T
    Tsim = temps(idT);
    display(Tsim)
        
    % And then, iterate over n_sims
    for idsim = 1:n_sims

        % Simulate AR(1) temperature
        tempar = zeros(timemax, 1);
        eta = randn(timemax, 1); % noise
        c = Tsim*(1 - estar1);
        tempar(1) = Tsim;
        % AR(1) temperature tie series
        for idx = 2:numel(tempar)
            tempar(idx) = c + estar1*tempar(idx-1) + sqrt(estvar)*eta(idx);
        end
        % Store temperature
        TempsAR{idT, idsim} = tempar;
        
        % Prepare forward Euler
        Rs = zeros(timemax, 1);
        Rval    = Rstart; % initial conditions
        time    = -dt;
        ii      = 0;
        
        % Do forward euler
        while time < timemax-1
            time    = time + dt;
            ii      = ii + 1;
            [dRdt, dVdt] = oerlemansModel(time, Rval, tempar(ii), par);
            Rval = Rval + dt*dRdt;
            % Store
            Rs(ii) = Rval;
        end % End forward Euler while loop
        % Store
        VresAR{idT, idsim} = OceanSurf*Vr(Rs);
    end % End idsim loop
end




%% Compute mean volume of fluctuating temperature simulations
meansAR = zeros(n_temps,1);
% Iterating over temperatures
for idT = 1:n_temps
        m = nan(n_sims, 1);
        % Iterating over each noisy simulation
        for idsim = 1:n_sims
            y = VresAR{idT, idsim};
            m(idsim) = mean( y((steps-last):end) );
        end
        % Store meansAR
        meansAR(idT) = mean(m);
end




%% Calculate the dV/dt matrix f, dV/dt = f(T, V)
% Oerlemans model with parameters fixed
ice = @(tx, Rx, Tx) oerlemansModel(tx, Rx, Tx, par);

% We define the interesting T and V region
% OceanSurf converts from volume to meters Sea Level Equivalent (m SLE)
Tvec = linspace(-4, 4, npoints);
Vvec = linspace(0, 8/OceanSurf,npoints);
% O/03 accepts R (radius) - we will convert
Rvec = nan(size(Vvec));

% Starting point for optimization; little magic
R0 = 100;
fminopt = optimset('Display', 'off');

% Find radius from volume
for ii = 1:numel(Vvec)
    af  = @(z) abs(Vvec(ii) - Vr(z));
    Rvec(ii) = fminsearch(af, R0, fminopt);
end

% Create meshgrids - need T, V and Rmat
[~, V] = meshgrid(Tvec, Vvec);
[T, Rmat] = meshgrid(Tvec, Rvec);

% Calculate f = dVdt
display('Calculating dV/dt')
f = nan(size(Rmat));
% We do this the old fashioned way - ice() is not vectorized
for ii = 1:size(Rmat, 1)
    for jj = 1:size(Rmat, 2)
        [dRdt, dVdt] = ice(0, Rmat(ii, jj), T(ii, jj));    
        f(ii, jj) = dVdt;
    end
end

% Convert to mm SLE 
f   = 1000*OceanSurf*f;
% Convert to m SLE
V   = V*OceanSurf;

% Temperature and Volume Step
dT = T(1,2) - T(1,1);
dV = V(2,1) - V(1,1);

% ------------ Aslaks figure ---------------------------------------------
fig002 = figure(002); hold on; box on; figset(fig002)
ph = pcolor(T, V, f);                  
set(ph, 'edgecolor', 'none');

% Use Aslaks colormap
cmap = hslcolormap(120, [0 nan 4]/6, 1, [0.1 .99 nan .99 0.1]);
colormap(cmap) 
clim = caxis;
caxis([-0.8 0.8]*max(abs(clim)));
hcb = colorbar;

% Contour of steady state
% We will use this handle later to draw a legend.
contour(T, V, f,[0 0],'k');


% Labels on axis and colorbar
xl = xlabel('Temperature $\bar{T}$ [$^{\circ}$C]'); textset(xl)
yl = ylabel('Ice Sheet Volume [m SLE]'); textset(yl)
%tt = title('$f = dV/dt$ [mm SLE yr$^{-1}$]'); textset(tt)
xlim([-2 2])
% Pretty colorbar
astring = 'dV/dt [mm SLE yr^{-1}]';
hcb.Label.String = astring;
set(hcb, 'FontSize', 14, 'FontName', 'Times New Roman')



%% ------------ Compute first, second and mixed derivatives ---------------
%
% https://se.mathworks.com/help/matlab/ref/gradient.html
% 
% [FX,FY,FZ,...] = gradient(F), where F has N dimensions, returns the N 
% components of the gradient of F. There are two ways to control the 
% spacing between values in F:
% 
%   1) A single spacing value, h, specifies the spacing between points in 
%   every direction.
%   2) N spacing values (h1,h2,...) specifies the spacing for each 
%   dimension of F. Scalar spacing parameters specify a constant spacing 
%   for each dimension. Vector parameters specify the coordinates of the 
%   values along corresponding dimensions of F. In this case, the length of 
%   the vector must match the size of the corresponding dimension.
% 
% *Note* The first output FX is always the gradient along the 2nd dimension 
% of F, going across columns. The second output FY is always the gradient 
% along the 1st dimension of F, going across rows. For the third output FZ 
% and the outputs that follow, the Nth output is the gradient along the Nth 
% dimension of F.

display('Doing derivatives')
% First derivative
[fT, fV] = gradient(f, dT, dV);
% Second and mixed
[fTT, fVT] = gradient(fT, dT, dV);
[fTV, fVV] = gradient(fV, dT, dV);




% Differentiate manually.

%    diff(X,N,DIM) is the Nth difference function along dimension DIM. 
%       If N >= size(X,DIM), diff returns an empty array.
 

% % df/dT
% fT = diff(f, 1, 2);
% fT = fT / dT;
% % d2f/dT
% fTT = diff(f, 2, 2);
% fTT = fTT / dT^2;
% % df/dV
% fV = diff(f, 1, 1);
% fV = fV / dV;
% % d2f/dV2
% fVV = diff(f, 2, 1);
% fVV = fVV / dV^2;
% 
% % d2f/dVdT
% fTV = diff(fT, 1, 1);
% fTV = fTV / dV;
% % d2f/dTdV
% fVT = diff(fV, 1, 2);
% fVT = fVT / dT;
% 
% 
% % Match dimensions
% fT = fT(3:end, 2:end);
% fTT = fTT(3:end, :);
% fV = fV(2:end, 3:end);
% fVV = fVV(:, 3:end);
% fTV = fTV(2:end, 2:end);
% fVT = fVT(2:end, 2:end);
% 
% % Matching f, volume and temperature matrices
% fm   = f(3:end, 3:end);
% T2m = T(3:end, 3:end);
% V2m = V(3:end, 3:end);


assert(max(max(abs(fVT - fTV))) < 1e-10, 'You must obey Schwarz theorem.')


%% Plot
% %------------- The PD wrt Temperature -----------------------------------
% fig003 = figure(003); box on; hold on; figset(fig003)
% ph = pcolor(T, V, fT);
% set(ph, 'edgecolor', 'none')
% colorbar;
% colormap(cmap) 
% clim = caxis;
% caxis([-0.4 0.4]*max(abs(clim)));
% xl = xlabel('Temperature $\bar{T}$ [$^{\circ}$C]'); textset(xl)
% yl = ylabel('Ice Sheet Volume [m SLE]'); textset(yl)
% xlim([-2 2])
% unit = '[mm yr $^{-1}$ $^{\circ}$C$^{-1}$]';
% tt = title(['$\partial f/ \partial T$', unit]); textset(tt)
% contour(T, V, f,[0 0],'k');


% % ------------- The PD wrt Volume -----------------------------------------
% fig004 = figure(004); box on; hold on; figset(fig004)
% ph = pcolor(T, V, fV);
% set(ph, 'edgecolor', 'none')
% colorbar;
% colormap(cmap) 
% clim = caxis;
% caxis([-0.4 0.4]*max(abs(clim)));
% xl = xlabel('Temperature $\bar{T}$ [$^{\circ}$C]'); textset(xl)
% yl = ylabel('Ice Sheet Volume [m SLE]'); textset(yl)
% xlim([-2 2])
% unit = ' [mm yr $^{-1}$ m $^{-1}$]';
% tt = title(['$\partial  f/ \partial  V$', unit]); textset(tt)
% contour(T, V, f,[0 0],'k');


% -------------- The second PD wrt Temperature twice ----------------------
fig005 = figure(005); box on; hold on; figset(fig005)
ph = pcolor(T, V, fTT);
set(ph, 'edgecolor', 'none')
colorbar;
colormap(cmap) 
clim = caxis;
caxis([-0.4 0.4]*max(abs(clim)));
xl = xlabel('Temperature $\bar{T}$ [$^{\circ}$C]'); textset(xl)
yl = ylabel('Ice Sheet Volume [m SLE]'); textset(yl)
xlim([-2 2])
unit = ' [mm yr $^{-1}$ $^{\circ}$C $^{-2}$]';
tt = title(['$ \partial^2 f / \partial T^2$', unit]); textset(tt)
contour(T, V, f,[0 0],'k');


% -------------- The second PD wrt Volume twice ---------------------------
fig006 = figure(006); box on; hold on; figset(fig006)
ph = pcolor(T, V, fVV);
set(ph, 'edgecolor', 'none')
colorbar;
colormap(cmap) 
clim = caxis;
caxis([-0.5 0.5]*max(abs(clim)));
xl = xlabel('Temperature $\bar{T}$ [$^{\circ}$C]'); textset(xl)
yl = ylabel('Ice Sheet Volume [m SLE]'); textset(yl)
xlim([-2 2])
unit = ' [mm yr $^{-1}$ m$^{-2}$]';
tt = title(['$\partial ^2 f /  \partial V^2$', unit]); textset(tt)
contour(T, V, f,[0 0],'k');


% --------------- The second PD wrt Temp first, then Volume ---------------
fig007 = figure(007); box on; hold on; figset(fig007)
ph = pcolor(T, V, fTV);
set(ph, 'edgecolor', 'none')
colorbar;
colormap(cmap) 
clim = caxis;
caxis([-0.5 0.5]*max(abs(clim)));
xl = xlabel('Temperature $\bar{T}$ [$^{\circ}$C]'); textset(xl)
yl = ylabel('Ice Sheet Volume [m SLE]'); textset(yl)
xlim([-2 2])
unit = ' [mm yr $^{-1}$ $^{\circ}$C $^{-1}$ m$^{-1}$]';
tt = title(['$d^2 f/ \partial  V  \partial T$', unit]); textset(tt)
contour(T, V, f,[0 0],'k');



%% Scatter
figure(002);
hold on,
scatter(temps, meanV);
scatter(temps, meansAR);

contour(T, V, f+fTT*arvar/2, [0,0], 'color', yellow, 'linewidth', 2)



%% Compute [(V - V*)^2](t) and [(V - V*)(T - T*)](t)

vv_fot  = cell(n_temps, n_sims);
vt_fot = cell(n_temps, n_sims);

for idT = 1:n_temps
    for idsim = 1:n_sims
        % Get temperature and volume af AR simulations
        Vhere = VresAR{idT, idsim};
        There = TempsAR{idT, idsim};
        % NORMALIZED
        %fprintf('***\n')
        %display(meansAR(idT))
        %display(temps(idT))
        %display( mean(Vhere((steps-last):end)))
        %display(mean(There))
        Vhere = Vhere - meansAR(idT);
        There = There - temps(idT);
        % Compute and store the values
        vv_fot{idT, idsim} = Vhere.^2;
        vt_fot{idT, idsim} = Vhere .* There;
        
    end
end

%% Plot!
cls = lines(n_temps);
close all;

% Legend strings
legstr = cell(n_temps,1);
plots = cell(n_temps,1);

for idx = 1:numel(temps)
    str = ['$\bar{T} = ', sprintf('%.1f', temps(idx)), '$'];
    legstr{idx} = str;
end


% Variance as a function of time
fig117 = figure(117); hold on; box on; figset(fig117)
subplot(211); hold on; box on
for idT = 1:n_temps
    for idsim = 1:n_sims
        pl = plot(vv_fot{idT, idsim}, 'Color', cls(idT, :));
        % Add this to plots
        if idsim == 1
            plots{idT, 1} = pl;
        end
    end
end



% Pretty
xl = xlabel('Time $t$ [years]'); textset(xl)
yl = ylabel('$[(V - \overline{V})^2](t)$'); textset(yl)
tt = title('$(V - \overline{V})^2$ as a function of time'); 
textset(tt)
ylim([-1 20])
xlim([0 2e4])

legvar = legend([plots{:}], legstr); legset(legvar)

% Covariance as a function of time
%fig119 = figure(119); hold on; box on; figset(fig119)
subplot(212); hold on; box on
for idT = 1:n_temps
    for idsim = 1:n_sims
        plot(vt_fot{idT, idsim}, 'Color', cls(idT, :))
    end
end

xl = xlabel('Time $t$ [years]'); textset(xl)
yl = ylabel('$[(V - \overline{V}) (T - \overline{T})](t)$'); textset(yl)
tt = title('$(V - \overline{V})(T - \overline{T})$ as a function of time'); 
textset(tt)
ylim([-20 20])
xlim([0 2e4])


%% Get values of fVV and fTV at the steady states

% Indices
tinds = zeros(numel(temps), 1);
for idx = 1:numel(temps)
    [~, index] = min(abs(Tvec - temps(idx)));
    tinds(idx) = index;
end

vinds = zeros(numel(meansAR), 1);
% Remember, we multiplied by OceanSurf
Vvec_alt = V(:, 1);
for idx = 1:numel(meansAR)
    [~, index] = min(abs(Vvec_alt - meansAR(idx)));
    vinds(idx) = index;
end


% Test
figure(002); hold on
for idx = 1:numel(temps)
    there = Tvec(tinds(idx));
    vhere = Vvec_alt(vinds(idx));
    
    scatter(there, vhere, 'kx')
end


%% Now, create vectors containing the values of fTT, fVV and fVT
fmat = zeros(numel(temps), 1);
fTTmat = zeros(numel(temps), 1);
fVVmat = zeros(numel(temps), 1);
fTVmat = zeros(numel(temps), 1);

for idx = 1:numel(temps)
    tind = tinds(idx);
    vind = vinds(idx);
    
    fmat(idx) = f(vind, tind);
    fTTmat(idx) = fTT(vind, tind);
    fVVmat(idx) = fVV(vind, tind);
    fTVmat(idx) = fTV(vind, tind);
end




%% Attempt to plot the above in a meaningful way

% Time axis
tax = (timemax-last):1:timemax;
ntimes = numel(tax);

% Save the plots
plotcells = cell(n_temps,1);


for idT=1:n_temps

    % Cell entry for this plot
    plotcells{idT,1} = figure(700+idT); hold on; box on; 
    figset(plotcells{idT, 1});

    % Plot f0 and fTT0
    subplot(2,1,1); hold on; box on
    ylim([-0.15 0.15])
    pf0_up = plot(tax, repmat(fmat(idT),  ntimes, 1), ...
        'linestyle', ':', 'linewidth', 2, 'color', blue);
    pfTT_up = plot(tax, repmat(fTTmat(idT)*arvar/2,ntimes, 1), ...
        'linestyle', ':', 'linewidth', 2, 'color', red);
    
    % Plot the same again...
    subplot(2,1,2); hold on; box on
    ylim([-0.15 0.15])
    pf0_dow = plot(tax, repmat(fmat(idT),  ntimes, 1), ...
        'linestyle', ':', 'linewidth', 2, 'color', blue);
    pfTT_dow = plot(tax, repmat(fTTmat(idT)*arvar/2,ntimes, 1), ...
        'linestyle', ':', 'linewidth', 2, 'color', red);

    % Compute means
    vv_means = zeros(1, n_sims);
    vt_means = zeros(1, n_sims);

    for idsim = 1:n_sims
        % (V - V0)^2 as a function of time
        vv = vv_fot{idT, idsim}(last:end);
        vv_means(idsim) = mean(vv);        
        % (T - T0)*(V - V0) as a function of time
        vt = vt_fot{idT, idsim}(last:end);
        vt_means(idsim) = mean(vt);
        
        % Plot
        subplot(2,1,1)
        plot(tax, 0.5*vv*fVVmat(idT));
        subplot(2,1,2)
        plot(tax, vt*fTVmat(idT));
    end
    
    % Compute means
    vv_mean = mean(vv_means);
    vt_mean = mean(vt_means);
    
    % Plot means times relevant derivative value
    subplot(2,1,1)
    phvv = plot(tax, repmat(0.5*fVVmat(idT)*vv_mean, ntimes, 1), 'k');
    subplot(2,1,2)
    phvt = plot(tax, repmat(fTVmat(idT)*vt_mean, ntimes, 1), 'k');
    
    % Pretty
    subplot(2,1,1)
    yl_up = ylabel('$f = dV/dt$'); textset(yl_up)    
    tt_up = title(['$[1/2 \times f_{VV}^0 \times(V - V_0)^2](t)$ ', ...
        'for $\bar{T}=$ ', sprintf('%.1f', temps(idT))]); textset(tt_up)
    
    subplot(2,1,2)
    xl_dow = xlabel('Simulation years'); textset(xl_dow)
    yl_dow = ylabel('$f = dV/dt$'); textset(yl_dow)
    
    tt_dow = title(['$[f_{TV}^0 \times (T - T_0)(V - V_0)](t)$ ', ...
        'for $\bar{T}=$ ', sprintf('%.1f', temps(idT))]); textset(tt_dow)
    
    % Upper legend
    legstrs_up = {'$f^0$',  ...
                  '$f_{TT}^0 \times Var(T) / 2$', ...
                  '$(1/2) \times \langle Var(V) \rangle \times f_{VV}^0$'};
    phs_up = [pf0_up, pfTT_up, phvv];
    subplot(2,1,1)
    leg_up = legend(phs_up, legstrs_up, 'location', 'northwest'); 
    legset(leg_up);
    
    %Lower legend
    legstrs_dow = {'$f^0$',  ...
                  '$f_{TT}^0 \times Var(T) / 2$', ...
                  '$\langle Cov(T, V) \rangle \times f_{TV}^0$'};
    phs_dow = [pf0_dow, pfTT_dow, phvt];
    subplot(2,1,2)
    leg_dow = legend(phs_dow, legstrs_dow, 'location', 'northwest'); 
    legset(leg_dow);
    
    % Verbosity
    vv_fvv = 0.5*vv_mean*fVVmat(idT);
    vt_ftv = vt_mean*fTVmat(idT);
    f0 = fmat(idT);
    
    fprintf('\n')
    fprintf('Tbar: %.1f\n', temps(idT))
    fprintf('vv_fvv: %.3s\n', vv_fvv)
    fprintf('vt_ftv: %.3s\n', vt_ftv)
    
    fprintf('vv_fvv/f0 * 100: %.3f\n', vv_fvv/f0*100)
    fprintf('vt_ftv/f0 * 100: %.3f\n', vt_ftv/f0*100)
end




%% Compute variance and covariance using MATLAB buil-ins

varmat = zeros(n_temps, n_sims);
covmat = zeros(n_temps, n_sims);

for idT = 1:n_temps
    for idsim = 1:n_sims
        
        % Temperature and
        Vhere = VresAR{idT, idsim}(last:end);
        There = TempsAR{idT, idsim}(last:end);
        
        covres = cov(There, Vhere);
        varres = var(Vhere);
        
        % Learning MATLABS cov()
        assert(abs(varres-covres(end, end)) < 1e-10);
        assert(covres(1,end) == covres(end,1));
        
        covmat(idT, idsim) = covres(1,end);
        varmat(idT, idsim) = varres;
    end
end

% Average
meanvar = mean(varmat, 2);
meancov = mean(covmat, 2);


% Multiply by the relevant derivatives
fVVterm = 0.5*meanvar.*fVVmat;
fTVterm = meancov.*fTVmat;

% 
fprintf('Term involving fTV in percents of f0:\n')
display(fTVterm ./ fmat * 100);
fprintf('Term involving fVV in percents of f0:\n')
display(fVVterm ./ fmat * 100);




%% Save pdf?
if save_pdf
    disp('Saving pdfs...')
    export_fig(fig117, [pdfpath, 'VarianceFunctionTime.pdf'])
end

if save_png
    disp('Saving pngs...')
    export_fig(fig117, [pngpath, 'VarianceFunctionTime.png'])
end

% -------------------------------------------------------------------------
fprintf('Done.\n')
toc
% -------------------------------------------------------------------------