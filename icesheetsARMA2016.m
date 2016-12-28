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
savefigures = true;
%savefigures = false;
savedata = true;
%savedata = false;
% 'test' is quite quick (~20 seconds), 'production' around 5 minutes
%thisrun = 'test';
thisrun = 'production';
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
tic
% -------------------------------------------------------------------------



% ------------- Pretty figures --------------------------------------------
fs = 18;
figset = @(f)   set(f, 'Color', 'w');
textset = @(h)  set(h, 'Fontsize', fs, 'Interpreter', 'Latex');
legset = @(l)   set(l, 'FontSize', fs, 'Interpreter', 'Latex');
hcbset = @(h) set(h, 'FontSize', fs, 'FontName', 'Times New Roman');
set(0,'DefaultAxesFontSize',12)
% -------------------------------------------------------------------------

% ------------- Prettier colors -------------------------------------------
cls = lines(5);
blue   = cls(1, :);
red    = cls(2, :);
yellow = cls(3, :);
purple = cls(4, :);
green  = cls(5, :);
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
addpath(genpath('MATLAB-Utilities'));
% Current working directory
currentwd = pwd;
% Parent folder
parentf = fileparts(currentwd);
% Where figures and data is stored
figpath = [parentf, '/LaTeX[IceSheets][2016]/gfx/'];
texpath = [parentf, '/LaTeX[IceSheets][2016]/tex/'];
datapath = [parentf, '/Data[IceSheets][2016]/'];
figformat = '.pdf';
% -------------------------------------------------------------------------


% ------------- Save figures? ---------------------------------------------
if savefigures == true;
    display('*--------------> Saving figures <--------------*')
end
% -------------------------------------------------------------------------
% Preserve order in the plots!
if savefigures == true
    set(0,'DefaultFigureWindowStyle','normal')
else
    set(0,'DefaultFigureWindowStyle','docked')
end
% -------------------------------------------------------------------------



% ------------- Save data? ------------------------------------------------
if savedata == true;
    display('*--------------> Saving data <-----------------*')
end
% -------------------------------------------------------------------------



% ------------- Load parameters -------------------------------------------
run('param2016')
% -------------------------------------------------------------------------


% ------------- Parameters for integration --------------------------------
dt           = 1.0;         % 1.0 in article
timemax      = 1e5;         % 5e4 in article
npoints      = 200;         % 2400 in article
% -------------------------------------------------------------------------


% ------------- Other parameters ------------------------------------------
if isequal(thisrun, 'test')
    n_sims      = 2; % Number of simulations per temperature
    temps_show = linspace(-2, 1.5, 3)';
    % Temps for steady numerical simulations plot
    temps = [-1.8; -0.8; 0; 1.5];
elseif isequal(thisrun, 'production')
    n_sims      = 5; % Number of simulations per temperature
    temps_show = linspace(-2, 1.5, 50)';
    % Temps for steady numerical simulations plot
    temps = [-1.8; -0.8; 0; 1.5];
end


n_temps     = numel(temps); % This many temperatures for OU tryouts
n_temps_show = numel(temps_show);
% Steady state is assumed after this simulation year
steady      = 5e4;
% Average the last 'last' values
last        = floor((timemax - steady)/dt);
% -------------------------------------------------------------------------



% ------------- Variance of ARMA(1, 1) process ---------------------------------
armavar = @(sig, phi, theta) sig^2 * ...
    (1 + theta^2 + 2*phi*theta) / (1 - phi^2);
% -------------------------------------------------------------------------


% ------------- From greenlandSigma.m -------------------------------------
try
    arma11results = load([datapath, 'arma11results.mat']);
catch FE
    display('No file ar1results.mat - run greenlandTemperatureARMA2016.m')
    % Can't proceed without these results
    rethrow(FE)
end
%display(ar1results.estar1)
%display(ar1results.estvar)
estar1 = arma11results.estar1;
estma1 = arma11results.estma1;
estvar = arma11results.estvar;

sqrtvar = sqrt(estvar);

armavar = armavar(sqrtvar, estar1, estma1);
% -------------------------------------------------------------------------





%% Integrate O/03 for different steady temperatures
display('Integrating steady temperatures')
% Volume as function of radius; 
Vr = @(R) (8/15).*pi.*sqrt(mu).*R.^(5/2) - (1/3).*pi.*s.*R.^(3);

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
        [dRdt, dVdt] = oerlemansClassic(time, Rval, Tsim, par);
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

        %{
        % Simulate AR(1) temperature
        tempar = zeros(timemax, 1);
        eta = randn(timemax, 1); % noise
        c = Tsim*(1 - estar1);
        tempar(1) = Tsim;
        
        
        % AR(1) temperature tie series
        for idx = 2:numel(tempar)
            tempar(idx) = c + estar1*tempar(idx-1) + sqrt(estvar)*eta(idx);
        end
        %}
        
        % Simulate ARMA(1, 1) temperature
        tempar = zeros(timemax, 1);
        eta = sqrtvar*randn(timemax+1, 1); % noise
        
        
        % Manual ARMA(1, 1)
        for idx = 2:numel(tempar)
            tempar(idx) = ...
                estar1*tempar(idx-1) + ...    % AR part of model
                estma1*eta(idx) + ...       % MA part of model
                eta(idx+1);                 % White nise
        end
        
        % We add the mean
        tempar = tempar + Tsim;
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
            [dRdt, dVdt] = oerlemansClassic(time, Rval, tempar(ii), par);
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

%% Compute <V - V*>, <(V - V*)^2> and <(V - V*)(T - T*)>
vv   = zeros(n_temps, n_sims);
vv2  = zeros(n_temps, n_sims);
vvtt = zeros(n_temps, n_sims);

for idT = 1:n_temps
    for idsim = 1:n_sims
        % Get temperature and volume af AR simulations
        Vhere = VresAR{idT, idsim};
        Vhere = Vhere((steps-last):end);
        There = TempsAR{idT, idsim};
        There = There((steps-last):end);
        % Compute means in the steady state
        Vmean = mean(Vhere);
        Tmean = mean(There);
        % Compute and store the values
        vv(idT, idsim) = mean(Vhere - Vmean);
        vv2(idT, idsim) = mean((Vhere - Vmean).^2);
        vvtt(idT, idsim) = mean((Vhere - Vmean).*(There - Tmean));
        
    end
end

%display(mean(vv, 2));
%display(mean(vv2, 2));
%display(mean(vvtt, 2));

% Save to LaTeX table
columnlabels = {...
    '$\bar{T}$ [$^{\circ}$C]', ...
    '$\bar{V}$', ...
    '$\langle (V - \bar{V})^2\rangle$', ...
    '$\langle (V - \bar{V})(T - \bar{T}) \rangle$', ...
    };
%    '$\hat{\sigma}(\bar{V}) / \bar{V}$'};
rowlabels = cellstr(num2str(temps, '%.2f'));
% Format the data
vfluc = zeros(n_temps, 4);
vfluc(:, 1) = temps;
vfluc(:, 2) = meansAR;
vfluc(:, 3) = mean(vv2, 2);
vfluc(:, 4) = mean(vvtt, 2);
%savedata(:, 5) = sqrt(mean(vv2, 2)) ./ meansAR;

if savedata == true
    savename = [texpath, 'fluctuations-ARMA.tex'];
    matrix2latex(vfluc, savename, ...    
    'columnLabels', columnlabels, ...
    'alignment', 'c', 'format', ...
    '%.2e', ...
    'size', 'large'); 
end
%'rowLabels', rowlabels, ...

% -------------------------------------------------------------------------
%% Repeat the above process for temps_show
Rresults_show = nan(numel(timeax), n_temps_show);
% Run the forward Euler on O/03
for idT = 1:n_temps_show
    Tsim = temps_show(idT);
    Rval    = Rstart; % initial conditions
    time    = -dt;
    ii      = 0;
    % Do forward euler
    while time < timemax
       time    = time + dt;
       ii      = ii + 1;
       [dRdt, dVdt] = oerlemansClassic(time, Rval, Tsim, par);
       Rval = Rval + dt*dRdt;
       % Store
       Rresults_show(ii, idT) = Rval;
     end % end while
end

% O/03 works in radius - we convert
Vresults_show = OceanSurf*Vr(Rresults_show);
% Compute mean volume of steady temperature simulations
meanV_show = zeros(n_temps_show, 1);
for idx = 1:n_temps_show
    meanV_show(idx) = mean(Vresults_show((steps-last):end, idx));
end

% Run the forward Euler on O/03 - with OU noise
VresAR_show = cell(n_temps_show, n_sims);
TempsAR_show = cell(n_temps_show, n_sims);
for idT = 1:n_temps_show % First, iterate over T
    Tsim = temps_show(idT);
    display(Tsim)
        
    % And then, iterate over n_sims
    for idsim = 1:n_sims

        %{
        % Simulate AR(1) temperature
        tempar = zeros(timemax, 1);
        eta = randn(timemax, 1); % noise
        c = Tsim*(1 - estar1);
        tempar(1) = Tsim;
        % AR(1) temperature tie series
        for idx = 2:numel(tempar)
            tempar(idx) = c + estar1*tempar(idx-1) + sqrt(estvar)*eta(idx);
        end
        %}
        
        % Simulate ARMA(1, 1) temperature
        tempar = zeros(timemax, 1);
        eta = sqrtvar*randn(timemax+1, 1); % noise
        
        
        % Manual ARMA(1, 1)
        for idx = 2:numel(tempar)
            tempar(idx) = ...
                estar1*tempar(idx-1) + ...    % AR part of model
                estma1*eta(idx) + ...       % MA part of model
                eta(idx+1);                 % White nise
        end
        
        % We add the mean
        tempar = tempar + Tsim;
        
        % Store temperature
        TempsAR_show{idT, idsim} = tempar;
        
        % Prepare forward Euler
        Rs = zeros(timemax, 1);
        Rval    = Rstart; % initial conditions
        time    = -dt;
        ii      = 0;
        
        % Do forward euler
        while time < timemax-1
            time    = time + dt;
            ii      = ii + 1;
            [dRdt, dVdt] = oerlemansClassic(time, Rval, tempar(ii), par);
            Rval = Rval + dt*dRdt;
            % Store
            Rs(ii) = Rval;
        end % End forward Euler while loop
        % Store
        VresAR_show{idT, idsim} = OceanSurf*Vr(Rs);
    end % End idsim loop
end

meansAR_show = zeros(n_temps_show,1);
% Iterating over temperatures
for idT = 1:n_temps_show
        m = nan(n_sims, 1);
        % Iterating over each noisy simulation
        for idsim = 1:n_sims
            y = VresAR_show{idT, idsim};
            m(idsim) = mean( y((steps-last):end) );
        end
        % Store meansAR
        meansAR_show(idT) = mean(m);
end
% -------------------------------------------------------------------------






%% Plot both fluctuating and steady temperatures as a function of time
fig001 = figure(001); hold on; box on; figset(fig001)
cls = lines(n_temps);

% First the steady...
phs = cell(n_temps, 1);
for idT = 1:n_temps
    plot(timeax, Vresults(:, idT), 'col', cls(idT,:), 'LineWidth', 2)
    %phs{idT} = ph;
end

% And now the fluctuating temperatures
for idsim = 1:max(2, n_sims)
    for idT = 1:n_temps
        plot(VresAR{idT, idsim}, 'col', cls(idT,:));
    end
end

% Pretty plot
%tt = title('Steady (thick line) and fluctuating temperature'); textset(tt)
xl = xlabel('Time [years]'); textset(xl)
yl = ylabel('Ice Sheet Volume [m SLE]'); textset(yl)

% Legend strings
legstr = cell(4,1);
for idx = 1:numel(temps)
    str = ['$\bar{V}(\bar{T} = ', sprintf('%.1f', temps(idx)), '$)'];
    legstr{idx} = str;
end

%legstr = [repmat('$\bar{T} = $', n_temps, 1), num2str(temps)];
lh = legend(legstr, 'location', 'northeast'); legset(lh)
ylim([0 7])


%% Calculate the dV/dt matrix f, dV/dt = f(T, V)
% Oerlemans model with parameters fixed
ice = @(tx, Rx, Tx) oerlemansClassic(tx, Rx, Tx, par);

% We define the interesting T and V region
% OceanSurf converts from volume to meters Sea Level Equivalent (m SLE)
Tvec = linspace(-4, 4, npoints);
Vvec = linspace(0, 8/OceanSurf,npoints);
% O/03 accepts R (radius) - we will convert
Rvec = nan(size(Vvec));

% Calculate R from V by optimisation - Equation (9)
% Volume as function of radius; 
Vr = @(R) (8/15).*pi.*sqrt(mu).*R.^(5/2) - (1/3).*pi.*s.*R.^(3);
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

% Contour folded with temperature fluctuation
%folding = normpdf(-14:dT:14, 0, arvar); 

%contour(T, V, imfilter(f, folding/sum(folding), ...
%    'replicate'), [0 0], 'LineColor', 'g', 'LineWidth', 2);

% Labels on axis and colorbar
xl = xlabel('Temperature $\bar{T}$ [$^{\circ}$C]'); textset(xl)
yl = ylabel('Ice Sheet Volume [m SLE]'); textset(yl)
%tt = title('$f = dV/dt$ [mm SLE yr$^{-1}$]'); textset(tt)
xlim([-2 2])
% Pretty colorbar
astring = 'dV/dt [mm SLE yr^{-1}]';
hcb.Label.String = astring;
set(hcb, 'FontSize', 14, 'FontName', 'Times New Roman')





%% Scatter steady state volume of constant and fluctuating T simulations
% fig003 = copy(fig002); hold on; box on; figset(fig003)
% Copy didn't work!!!!?!
fig003 = figure(003); hold on; box on; figset(fig003)

% Redo figure 2
ph = pcolor(T, V, f);                  
colormap(cmap) 
clim = caxis;
caxis([-0.8 0.8]*max(abs(clim)));
set(ph, 'edgecolor', 'none');
hcb = colorbar;

% Contour of steady state
[C0, c0cont] = contour(T, V, f,[0 0],'k');

% Scatter 
s1 = scatter(temps_show, meanV_show, ...
    'markeredgecolor', green, 'markerfacecolor', 'none');
s2 = scatter(temps_show, meansAR_show, ...
    'markeredgecolor', blue, 'markerfacecolor', blue);

% Legend
shs = [c0cont, s1, s2];
legstr = {'$f=0$', ...
    'Constant $\bar{T}$', ...
    'Fluctuating $T_t$'};
leg3 = legend(shs, legstr); legset(leg3)

% Pretty plot
xl = xlabel('Temperature $\bar{T}$ [$^{\circ}$C]'); textset(xl)
yl = ylabel('Ice Sheet Volume [m SLE]'); textset(yl)
%tt = title('$f := dV/dt$ and results of numerical simulations (dots)'); textset(tt)
xlim([-2 2])
% Pretty colorbar
astring = 'dV/dt [mm SLE yr^{-1}]';
hcb.Label.String = astring;
set(hcb, 'FontSize', 14, 'FontName', 'Times New Roman')





%% Compute ACF(V) for ALL volume simulations
acfall = cell(n_temps, n_sims);
cctall = cell(n_temps, n_sims);

% nlags lags - ylast years
nlags = floor(last/2);


for idT = 1:n_temps
    for idsim = 1:n_sims
        vh = VresAR{idT, idsim}(end-last+1:end);
        [acf, lags, ~] = autocorr(vh, nlags);
        %plot(lags, acf)
        acfall{idT, idsim} = acf;
        cctall{idT, idsim} = find(acf < exp(-1), 1);
    end
end

%% Calculate correlation time of the ice sheet
fig004 = figure(004); hold on; box on; figset(fig004)

colors = lines(n_temps);
idc = 0;

phs = cell(n_temps, 1);

for idT = 1:n_temps
    idc = idc + 1;
    for idsim = 1:n_sims
        acf = acfall{idT, idsim};
        ph = plot(lags, acf, 'color', colors(idc,:));
        % Append to cell of legend handles
        if idsim == 1
            phs{idT} = ph;
        end
    end
end


% Plot 1/e
oe = repmat(exp(-1), nlags+1, 1);
poe = plot(lags, oe, 'color', 0.4*[1 1 1], 'LineWidth', 1.7);

%tt = title('ACF(k) of steady state volume with fluctuating temperature'); 
% textset(tt)
xl = xlabel('$k$ (Lag)'); textset(xl)
yl = ylabel('ACF(k)'); textset(yl)

% Legend strings... I miss Python sometimes.
legstrs = cell(n_temps+1, 1);
legstrs{1} = '$1/e$';
for idx = 2:n_temps+1
    legstrs{idx} = ['$\bar{T} =' num2str(temps(idx-1)), '$'];
end
% The legend itself
l1 = legend([poe phs{:}], legstrs); legset(l1)

% Write mean, max, min, std of cctall
cctmat = cell2mat(cctall);







%% Scatter the fluctuating SIMULATIONS on the figure
fig006 = copy(fig002); figset(fig006)
% Only scatter lastsc points - there are many!
lastsc = 1000;

for idT = 1:n_temps
    for idsim = 1:n_sims
        % Pick values from AR simulations
        t = TempsAR{idT, idsim};
        v = VresAR{idT, idsim};
        % Only keep the last from steady state
        t = t((steps-lastsc):end);
        v = v((steps-lastsc):end);
        scatter(t, v, 'k.')
    end
end


%tt = title('$V_t$ scattered vs. $T_t$ for different realisations');
%textset(tt)

s1 = scatter(temps, meanV, ...
    'markeredgecolor', green, 'markerfacecolor', green);
s2 = scatter(temps, meansAR, ...
    'markeredgecolor', blue, 'markerfacecolor', blue);

xlim([-2 2])
ylim([0 8])


%% OK - draw probability distributions instead
fig007 = copy(fig002); figset(fig007); hold on; box on

cls = lines(3);

lastsc=5e4;

for idT = 1:n_temps
idsim=randi(n_sims);
    % Pick values from AR simulations
    t = TempsAR{idT, idsim};
    v = VresAR{idT, idsim};
    % Only keep the last from steady state
    t = t((steps-lastsc):end);
    v = v((steps-lastsc):end);

    % Now we do the normpdf
    % Here, the x- and y-axes on which to plot
    there = mean(t);
    vhere = mean(v);
    tax = linspace(there-2, there+2, 501);
    % Uneven number of points ensure we can see the delta function
    vax = linspace(vhere-1, vhere+1.5, 101);

    tpdf = normpdf(tax, there, std(t));
    vpdf = normpdf(vax, vhere, std(v));

    plot(tax, tpdf+vhere, 'k')% 'color', cls(1,:))
    plot(tax, vhere*ones(size(tax)), 'k:')

    plot(vpdf/max(vpdf)+there, vax, 'k')%'color', cls(2, :))
end

scatter(temps, meansAR, ...
    'markeredgecolor', blue, 'markerfacecolor', blue);


xlim([-2 2])
ylim([0 8])

%tt = title('PDF of $V_t$ and $T_t$ for different realisations');
%textset(tt)


%% Find out how close <Vt> is to \bar{V}
%  Run the analysis in intervals of the mean correlation time
meancct = floor(mean(cctmat(:)));
stdcct  = floor(std(cctmat(:)));

%meancct = 1e4;

% This many intervals
intervals = floor(last/meancct);
confints = [0.7, 0.8, 0.9, 0.95, 0.99];

devmat = zeros(n_temps, numel(confints), n_sims, intervals);
widthmat = zeros(n_temps, numel(confints), n_sims, intervals);
devmat(:) = nan;

% Iterate over temps
for idT = 1:n_temps
    % Iterate over the confidence intervals
    for idconf = 1:numel(confints)
        confval = 1-confints(idconf);
        % Iterate over simulations
        for idsim = 1:n_sims
            % Pick values from AR simulations
            v = VresAR{idT, idsim};
            % Only keep the last from steady state
            v = v((steps-last):end);
            % Iterating over the various intervals
            for idint = 1:intervals
                % This subinterval will be interesting
                vhere = v((idint-1)*meancct+1:idint*meancct);
                vmean = mean(vhere);
            
                % Do the above with ksdensity()
                [fks, xks] = ksdensity(vhere);
                %plot(xks, fks, 'Color', colors(ii,:))

                % Construct credible intervals
                ct = cumtrapz(xks, fks);
                %display(ct)
                imin = find(ct < confval/2, 1, 'last');
                imax = find(ct > 1 - confval/2, 1, 'first') + 1;
                
                % The minimum and maximum deviation at the given confint
                minp = xks(imin);
                maxp = xks(imax);
                                
                % Distance
                d = maxp - minp;
                
                %fprintf('T: %.2f; Conf.: %.2f\n', temps(idT), confint)
                %fprintf('V: %.2f; d: %.2e; d/m: %.2e\n', ...
                %    vmean, d, d/vmean)
                %fprintf('imin: %i, imax: %i', imin, imax)
                %fprintf('\n')
                                
                % Width of interval divided by mean
                deviation = d/vmean;
                
                % Store
                devmat(idT, idconf, idsim, idint) = deviation;
                widthmat(idT, idconf, idsim, idint) = d;
            end
        end
    end
end


%% Plot the deviations
fig008 = figure(008); hold on; box on; figset(fig008)
cls = lines(n_temps);

dev_mean_intervals = mean(devmat, 4);
dev_mean_sims = mean(dev_mean_intervals, 3);

width_mean_intervals = mean(widthmat, 4);
width_mean_sims = mean(width_mean_intervals, 3);


dev_std_sims = zeros(size(dev_mean_sims));
% Standard deviations. For loop is deemed safer for the 4D array...
for idT=1:n_temps
    for idconf=1:numel(confints)
    v = devmat(idT, idconf, :, :);
    dev_std_sims(idT, idconf) = std(v(:));
    end
end

stderrs = dev_std_sims/sqrt(intervals*n_sims);

% Check to see that the above was done correctly
%assert(isequal(size(deviations), [n_temps numel(confints)]));

phs = cell(n_temps, 1);
for idT = 1:n_temps
    devs = dev_mean_sims(idT, :);
    stderr = stderrs(idT, :);
    %
    p = errorbar(confints, 100*devs, 100*stderr, 100*stderr, ...
        'col', cls(idT,:), 'LineWidth', 1.3);
    phs{idT} = p;
end


% Pretty plot
%tt = title('Magnitude of Steady State Fluctuations'); textset(tt)
xl = xlabel('Credibility Level of SSIW'); textset(xl)
yl = ylabel(['$\langle$ SSIW$^{(i)}$ / $\langle V_t^{(i)}\rangle \rangle \times 100\%$']); textset(yl)
legstr = [repmat('$\bar{T} = $', n_temps, 1), num2str(temps)];
lh = legend([phs{:}], legstr, 'location', 'northwest'); legset(lh)

%% Construct table-ifyable matrix whith these values:
%
% ssiwmat: 
% Temp  CL1         CL2         ...     CLN
%-------------------------------------------------------------------------
% T1    Width(1,1)  Width(1,2)  ...     Width(1,N)
% .     .
% .         .
% TN            .
%
% ssiwsemat:
% Temp  CL1         CL2         ...     CLN
%-------------------------------------------------------------------------
% T1    +/-         +/1
% .     .
% .         .
% TN            .
ssiwmat = nan(n_temps+1, 1+numel(confints));
ssiwmat(1, 2:end) = confints';
ssiwsemat = ssiwmat;

for idx = 2:n_temps+1
    devs = dev_mean_sims(idx-1, :);
    stderr = stderrs(idx-1, :);
    
    % Interval widths - first, temperature
    ssiwmat(idx, 1) = temps(idx-1);
    ssiwmat(idx, 2:end) = 100*devs;
    
    % Interval widths - first, temperature
    ssiwsemat(idx, 1) = temps(idx-1);
    ssiwsemat(idx, 2:end) = 100*stderr;
end


columnlabels = cell(numel(confints)+1, 1);
columnlabels{1} = '$\bar{T}$';

for idx = 1:numel(confints)
    columnlabels{idx+1} = num2str(confints(idx));
end


if savedata == true
    savename1 = [texpath, 'ssiw-ARMA.tex'];
    matrix2latex(ssiwmat(2:end, :), savename1, ...    
    'columnLabels', columnlabels, ...
    'alignment', 'c', 'format', ...
    '%.2f', ...
    'size', 'large'); 

    savename2 = [texpath, 'ssiwse-ARMA.tex'];
    matrix2latex(ssiwsemat(2:end, :), savename2, ...    
    'columnLabels', columnlabels, ...
    'alignment', 'c', 'format', ...
    '%.2f', ...
    'size', 'large'); 
end



%% Approximate the steady state of O/03 with fluctuating temperature

% Differentiate. We are now missing two rows on the left.
d2f = diff(f, 2, 2);
d2f = d2f / dT^2;
% Matching f, volume and temperature matrices
fm   = f(:, 3:end);
T2m = T(:, 3:end);
V2m = V(:, 3:end);


fig009 = figure(009); hold on; box on; figset(fig009)
ph = pcolor(T2m, V2m, fm);                  
set(ph, 'edgecolor', 'none');

% Use Aslaks colormap
cmap = hslcolormap(120, [0 nan 4]/6, 1, [0.1 .99 nan .99 0.1]);
colormap(cmap) 
clim = caxis;
caxis([-0.8 0.8]*max(abs(clim)));
hcb = colorbar;

% Pretty colorbar
astring = 'dV/dt [mm SLE yr^{-1}]';
hcb.Label.String = astring;
set(hcb, 'FontSize', 14, 'FontName', 'Times New Roman')


% Scatter 
s1 = scatter(temps_show, meanV_show, 'k+');
%    'markeredgecolor', 'k', 'markerfacecolor', 'none', 'linewidth', 1);

% Now, prepare to do the magic stuff
[C3, d2fcont] = contour(T2m, V2m, fm+d2f*armavar/2, ...
    [0 0], 'LineColor', yellow, 'LineWidth', 2);

% Also, scatter the fluctuating temperatures on the color map
sc = scatter(temps_show, meansAR_show, 'ko');
%    'markeredgecolor', yellow, 'markerfacecolor', 'none', 'linewidth', 2);


% Contour of steady state
[C1, odecont] = contour(T2m, V2m, fm,[0 0],'k');

% Contour folded with temperature fluctuation
%folding = normpdf(-14:dT:14, 0, sqrt(arvar)); 





% Legend
shs = [odecont, d2fcont, s1, sc];
legstr = {'$f=0$', ...
    '$f^* + \sigma_T^2/2 \times f^*_{TT}= 0$', ...
    'Constant $\bar{T}$', ...
    'Fluctuating $T_t$'};

leg3 = legend(shs, legstr); legset(leg3)




% Legend
%legstrs = {'$f^* + \sigma_T^2/2 \times f^*_{TT}= 0$', ...
%    'Numerical Simulations'};
%l1 = legend([d2fcont sc], legstrs); legset(l1)

xl = xlabel('Temperature $\bar{T}$ [$^{\circ}$C]'); textset(xl)
yl = ylabel('Ice Sheet Volume [m SLE]'); textset(yl)

xlim([-2 2])
ylim([0 8])





%% Illustrate the nonlinearity of dV/dt

% We plot for volumes in valwant
valwant = [2 3 4 5 6 7 ];
nvals   = numel(valwant);
indics  = zeros(size(valwant));

% Find the indices in V corresponding to {Vi | Vi in valwant}
for ii = 1:numel(valwant)
    [~, idx] = min(abs(V(:, 1) - valwant(ii)));
    indics(ii) = idx;
end

% 1) Pick the middle curve; for illustration purposes
middle = indics(ceil(numel(indics) / 2));

% Figure: dV/dt(T) and some illustration
fig010 = figure(010); hold on; box on; figset(fig010)
legs = zeros(nvals, 1);

% Plot dV/dt(T) for different volumes
for ii = 1:nvals
    fT = f(indics(ii),:);
    legs(ii) = V(indics(ii), 1);
    p = plot(Tvec, fT);
    % We want to emphasize the middle one later
    if indics(ii) == middle
        set(p, 'LineWidth', 2.5)
        keepcol = get(p, 'Color');
    end
end


% Pretty plot
xl = xlabel('Temperature $\bar{T}$ [$^{\circ}$C]'); textset(xl)
yl = ylabel('$\dot{V}$ [mm SLE yr$^-1$]'); textset(yl)
l1 = legend(num2str(legs, '%1.1f'), 'Location', 'SouthEast'); legset(l1)
legpos = l1.Position;
% Write 'Volume' over the legend box
text(2.75, -7.8, ...
    'Volume', 'Fontsize', 14, 'Interpreter', 'Latex',...
    'BackgroundColor', 'w');


% Now do the picture-in-picture. Values are chosen for visibility.
xpos = Tbar-1.2*sqrt(armavar);    % Left corner
w = 2.4*sqrt(armavar);            % Width
ypos = -2.5;                    % Lower edge
h = 3;                          % Height
% Draw rectangle around
rectangle('Position', [xpos ypos w h], 'LineStyle', '-.')

% We are plotting this curve
interestvec = f(middle,:);
% Make box for plot-in-plot
ax2 = axes('Position', [0.2, 0.15  0.45 0.45]);
hold on; box on
plot(ax2, Tvec, interestvec, 'Color', keepcol)
xlim([xpos, xpos+w])
ylim([ypos, ypos+h])


% Plot 0, +/- sigma
[~, ind0] = min(abs(Tvec - Tbar));
[~, indP] = min(abs(Tvec - Tbar - sqrt(armavar)));
[~, indM] = min(abs(Tvec - Tbar + sqrt(armavar)));

%hold on
scatter(Tvec(ind0), interestvec(ind0), 'ko')
scatter(Tvec(indP), interestvec(indP), 'ko')
scatter(Tvec(indM), interestvec(indM), 'ko')

% Calculate f(x) = ax + b;
% a = (y2 - y1) / (x2 - x1)
ap = (interestvec(indP) - interestvec(indM)) / (Tvec(indP) - Tvec(indM));
bp = interestvec(indP) - ap * Tvec(indP);
fp = @(x) ap.*x + bp;
plot(Tvec, fp(Tvec), 'k-.')

% %http://www.mathworks.com/matlabcentral/answers/4515-removing-ticks
set(gca, 'xtick', [])
set(gca, 'ytick', [])

% Arrows
% http://www.mathworks.com/matlabcentral/fileexchange/10656-data-space-to-figure-units-conversion

yl = ylim;
ydist = yl(2) - yl(1);
sigstr = sprintf('$\\sigma$');%, 'Fontsize', fs, 'Interpreter', 'Latex');

xa1 = [Tbar - sqrt(armavar); Tbar-0.005];
yp1 = (yl(1) + 0.6*(yl(2) - yl(1)));
ya1 = interestvec(ind0)*[1 1];

[xaf1, yaf1] = ds2nfu(xa1, ya1);
annotation('doublearrow', xaf1, yaf1)
text(Tbar - 0.5*sqrt(armavar), yp1-0.05*ydist, sigstr, 'Fontsize', 18, 'Interpreter', 'Latex');

xa2 = [Tbar+0.03; Tbar + sqrt(armavar)];
yp2 = yl(1) + 0.5*(yl(2) - yl(1));
ya2 = interestvec(ind0)*[1 1];

[xaf2, yaf2] = ds2nfu(xa2, ya2);
annotation('doublearrow', xaf2, yaf2 )
text(Tbar + 0.5*sqrt(armavar), yp2+0.05*ydist, sigstr, 'Fontsize', 18, 'Interpreter', 'Latex');


% Vertical dashed lines
l1 = line(Tvec(indM)*[1 1], [ya1(1) interestvec(indM)]);
set(l1, 'Color', 'k', 'LineStyle', '-.')
%
l2 = line(Tvec(indP)*[1 1], [ya1(2) interestvec(indP)]);
set(l2, 'Color', 'k', 'LineStyle', '-.')




%% Save figures?
if savefigures
    fprintf('Saving figures...\n');
    export_fig(fig001, [figpath, 'Simulations-ARMA-2016',    figformat])
    export_fig(fig003, [figpath, 'SteadyStates-ARMA-2016',   figformat])
    export_fig(fig004, [figpath, 'ACF(V)-ARMA-2016',         figformat])
    export_fig(fig006, [figpath, 'SimScatter-ARMA-2016',     figformat])
    export_fig(fig007, [figpath, 'SimProbability-ARMA-2016', figformat])
    export_fig(fig008, [figpath, 'SteadyInterval-ARMA-2016', figformat])
    export_fig(fig009, [figpath, 'SteadyApprox-ARMA-2016',   figformat])
    export_fig(fig010, [figpath, 'MassBalance-ARMA-2016',    figformat])
end


% -------------------------------------------------------------------------
fprintf('\n')
toc
% -------------------------------------------------------------------------



%% Print out data
fprintf('\n\n    ****************** Computed Values *****************\n\n')

fprintf('----------- ACF(volume) ---------------------------\n')
fprintf('Mean CCT: %.3f\n', mean(cctmat(:)))
fprintf('Min CCT: %.3f\n', min(cctmat(:)))
fprintf('Max CCT: %.3f\n', max(cctmat(:)))
fprintf('STD CCT: %.3f\n', std(cctmat(:)))

% Compare CCT(T) and CCT(V)
tauTmin = min(cctmat(:));

fprintf('\n')
%fprintf('CCT(T) rounded up to nearest year: %.2f\n', ceil(tauT))
fprintf('min([CCT(V)]) rounded up to nearest year: %.2f\n', ceil(tauTmin));
%fprintf('ceil(min[CCT(V)]) / ceil(CCT[T]): %.2f\n', ...
%    ceil(tauTmin)/ceil(tauT));

format bank % two decimal points; 'currency' format
fprintf('\n\n----------- SSIW ---------------------------\n')
fprintf('First row: confidence levels\n')
fprintf('First column: temperature\n\n')
fprintf('SSIW:\n')
display(ssiwmat)
fprintf('\n')
fprintf('Standard error of SSIW:\n')
display(ssiwsemat)
format short

format short % default