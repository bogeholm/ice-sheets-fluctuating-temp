%% Integrate Hans Oerlemans' (2003) models with noise
%
% Troels B. Mikkelsen - bogeholm@nbi.ku.dk
% 2015 - 2016


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
% Save data?
%savedata = true;
savedata = false;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 'test' is quite quick run (~20 seconds on dual core i7), 
% 'production' around 1 minute
%thisrun = 'test';
thisrun = 'production';
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Turn off warnings for ds2nfu
warning('off', 'MATLAB:nargchk:deprecated')
% -------------------------------------------------------------------------




% ------------- Verbosity ---------------------------------------------
if savedata == true
    disp('*--------------> Saving data <-----------------*')
end
% -------------------------------------------------------------------------


% ------------- Parameters for integration --------------------------------
dt           = 1.0;         % 1.0 in article
timemax      = 1e5;         % 1e5 in article
% -------------------------------------------------------------------------


% ------------- Other parameters ------------------------------------------
if isequal(thisrun, 'test')
    n_sims      = 2; % Number of simulations per temperature
    % Temps for steady numerical simulations plot
    temps = [-1.5; 0; 1.5; 3];
    temps_show = temps;
    npoints      = 300;         % 
elseif isequal(thisrun, 'production')
    n_sims      = 5; % Number of simulations per temperature
    % Temps for steady numerical simulations plot
    temps = [-1.5; 0; 1.5; 3];
    npoints      = 1200;         % 1200 in article
    temps_show = temps;
end


n_temps     = numel(temps); % This many temperatures for OU tryouts
n_temps_show = numel(temps_show);
% Steady state is assumed after this simulation year
steady      = 5e4;
% Average the last 'last' values
last        = floor((timemax - steady)/dt);
% -------------------------------------------------------------------------



% ------------- Variance of AR(1) process ---------------------------------
%https://en.wikipedia.org/wiki/Autoregressive_model#
% Example:_An_AR.281.29_process
%
% x(t+1) = a1*x(t) + s*eta => var(xt) = s^2/(1 - a1^2)
arvar_func = @(a1, s) s^2 / (1 - a1^2);
% -------------------------------------------------------------------------


% ------------- From greenlandSigma.m -------------------------------------
try
    ar1results = load([datapath, 'ar1results.mat']);
catch FE
    disp('No file ar1results.mat. Please run greenlandTemperature.m')
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
disp('Integrating steady temperatures')
% Volume as function of radius; 
Vr = @(R) iceTotal(R, par);


Rstart = 7e5;   % Initial volume; inspect results to ensure sanity
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
disp('Integrating OU fluctuating temperatures')

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


meansAR= zeros(n_temps, 1);
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
       [dRdt, dVdt] = oerlemansModel(time, Rval, Tsim, par);
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
            [dRdt, dVdt] = oerlemansModel(time, Rval, tempar(ii), par);
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
    str = ['$V_t(\bar{T} = ', sprintf('%.1f', temps(idx)), '$)'];
    legstr{idx} = str;
end

%legstr = [repmat('$\bar{T} = $', n_temps, 1), num2str(temps, 2)];
lh = legend(legstr, 'location', 'northwest'); legset(lh)
ylim([2 10])



%% Calculate the dV/dt matrix f, dV/dt = f(T, V)
% Oerlemans model with parameters fixed
ice = @(tx, Rx, Tx) oerlemansModel(tx, Rx, Tx, par);

% We define the interesting T and V region
% OceanSurf converts from volume to meters Sea Level Equivalent (m SLE)
Tvec = linspace(-2, 4, npoints);
Vvec = linspace(0, 11/OceanSurf,npoints);
% O/03 accepts R (radius) - we will convert
Rvec = nan(size(Vvec));

% Calculate R from V by optimisation - Equation (9)
% Volume as function of radius; 
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
disp('Calculating dV/dt')
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
%xlim([-4 6])
% Pretty colorbar
astring = 'dV/dt [mm SLE yr^{-1}]';
hcb.Label.String = astring;
set(hcb, 'FontSize', fs, 'FontName', 'Times New Roman')





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
set(hcb, 'FontSize', fs, 'FontName', 'Times New Roman')


% Scatter 
s1 = scatter(temps_show, meanV_show, 'k+');
%    'markeredgecolor', 'k', 'markerfacecolor', 'none', 'linewidth', 1);

% Now, prepare to do the magic stuff
[C3, d2fcont] = contour(T2m, V2m, fm+d2f*arvar/2, ...
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
    '$f^0 + \sigma_T^2/2 \times f^0_{TT}= 0$', ...
    'Constant $\bar{T}$', ...
    'Fluctuating $T_t$'};

leg3 = legend(shs, legstr); legset(leg3)




% Legend
%legstrs = {'$f^* + \sigma_T^2/2 \times f^*_{TT}= 0$', ...
%    'Numerical Simulations'};
%l1 = legend([d2fcont sc], legstrs); legset(l1)

xl = xlabel('Temperature $\bar{T}$ [$^{\circ}$C]'); textset(xl)
yl = ylabel('Ice Sheet Volume [m SLE]'); textset(yl)

xlim([-2 4])
ylim([1 11])



%% Now plot simulations and dV/dt on same plot
subplotfont = 10;
% We might be able to do this with copyobj, ... but do it brute force here.
fig011 = figure(011); figset(fig011)

% Black and grey color
black = [0 0 0];
grey = 0.6*[1 1 1];





spleft = subplot(121); hold on; %box on
% First the steady...
phs = cell(n_temps, 1);
for idT = 1:n_temps
    plot(timeax, Vresults(:, idT), 'col', black, 'LineWidth', 1)
    %phs{idT} = ph;
end

% tick size
set(gca, 'FontSize', subplotfont);



% And now the fluctuating temperatures
for idsim = 1:max(2, n_sims)
    for idT = 1:n_temps
        plot(VresAR{idT, idsim}, 'col', grey);
    end
end

% Pretty plot
%textset = @(h)  set(h, 'Fontsize', fs, 'Interpreter', 'Latex');
%tt = title('Steady (thick line) and fluctuating temperature'); textset(tt)
xl = xlabel('Time [years]'); %textset(xl)
set(xl, 'Fontsize', subplotfont, 'Interpreter', 'Latex');
yl = ylabel('Ice Sheet Volume [m SLE]'); %textset(yl)
set(yl, 'Fontsize', subplotfont, 'Interpreter', 'Latex');

% Legend strings
legstr = cell(4,1);
for idx = 1:numel(temps)
    str = ['$V_t(\bar{T} = ', sprintf('%.1f', temps(idx)), '$)'];
    legstr{idx} = str;
end

%legstr = [repmat('$\bar{T} = $', n_temps, 1), num2str(temps, 2)];
%lh = legend(legstr, 'location', 'northwest'); legset(lh)


% The other plot
spright = subplot(122); hold on; box on

ph = pcolor(T2m, V2m, fm);                  
set(ph, 'edgecolor', 'none');

% tick size
set(gca, 'FontSize', subplotfont);

colormap(cmap) 
clim = caxis;
caxis([-0.8 0.8]*max(abs(clim)));

hcb = colorbar;

% Pretty colorbar
astring = 'dV/dt [mm SLE yr^{-1}]';
hcb.Label.String = astring;
set(hcb, 'FontSize', subplotfont, 'FontName', 'Times New Roman')

% Scatter 
s1 = scatter(temps_show, meanV_show, 'k+');
%    'markeredgecolor', 'k', 'markerfacecolor', 'none', 'linewidth', 1);

% Now, prepare to do the magic stuff
[C3, d2fcont] = contour(T2m, V2m, fm+d2f*arvar/2, ...
    [0 0], 'LineColor', yellow, 'LineWidth', 2);

% Also, scatter the fluctuating temperatures on the color map
sc = scatter(temps_show, meansAR_show, 'ko');
%    'markeredgecolor', yellow, 'markerfacecolor', 'none', 'linewidth', 2);

% Contour of steady state
[C1, odecont] = contour(T2m, V2m, fm,[0 0],'k');

% Legend


xl = xlabel('Temperature $\bar{T}$ [$^{\circ}$C]'); % textset(xl)
set(xl, 'Fontsize', subplotfont, 'Interpreter', 'Latex');

% % Vertical lines
% subplot(122); hold on
% for idx = 1:numel(temps)
%     v = meanV(idx);
%     t = temps(idx);
%     line(t*[1 1], [1 v], 'Color', cls(idx, :), 'LineStyle', ':')
% end


% Adjust limits
subplot(121);
ylim([1 11]);

subplot(122);
ylim([1 11]);
set(gca, 'ytick', []);

% Make right subplot wider
% https://se.mathworks.com/matlabcentral/answers/...
% 119956-how-can-i-make-subplots-larger
shrink = 0.1; % Move spleft this much left, right this much right
enlarge = 0.8; % enlarge spright this much
flatten = 0.3; % flatten subplots this amount

posleft = get(spleft, 'Position'); % gives the position of current sub-plot
new_posleft = posleft +[0 0 -shrink -0.3];
set(spleft, 'Position', new_posleft ); % set new position of current sub - plot

posright = get(spright, 'Position'); % gives the position of current sub-plot
new_posright = posright +[-shrink*(1 + enlarge) 0 shrink -0.3];
set(spright, 'Position', new_posright ); % set new position of current sub - plot



% Normalized figure units for left axes - steady
[xlc, ylc] = ds2nfu(spleft, repmat(timemax, n_temps, 1), meanV);
% Right axes, same
[xrc, yrc] = ds2nfu(spright, temps, meanV);

% Draw these lines
for idx = 1:n_temps
    x = [xlc(idx) xrc(idx)];
    y = [ylc(idx), yrc(idx)];
    annotation('line', x, y, 'Color', black, 'Linewidth', 1, 'Linestyle', ':');
end


% Normalized figure units for left axes - noisy
flucstart = 0.5*timemax;
[xln, yln] = ds2nfu(spleft, repmat(flucstart, n_temps, 1), meansAR);
% Right axes, same
[xrn, yrn] = ds2nfu(spright, temps, meansAR);

% Draw these lines
for idx = 1:n_temps
    x = [xln(idx) xrn(idx)];
    y = [yln(idx), yrn(idx)];
    annotation('line', x, y, 'Color', grey, 'Linewidth', 1, 'Linestyle', ':');
end


% Legend
shs = [odecont, d2fcont, s1, sc];
legstr = {'$f=0$', ...
    '$f^0 + \frac{\sigma_T^2}{2} f^0_{TT}= 0$', ...
    'Constant $\bar{T}$', ...
    'Fluctuating $T_t$'};

varlegset = @(l)   set(l, 'FontSize', subplotfont, 'Interpreter', 'Latex');
leg011 = legend(shs, legstr); varlegset(leg011)



% tic
% desktop = '/Users/bogeholm/Desktop/';
% print(fig011, [desktop, 'Sim+Approx-2016.eps'], '-depsc')
% toc

%% Illustrate the nonlinearity of dV/dt

% We plot for volumes in valwant
valwant = [3; 4; 5; 6; 7; 8];
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
orgplot = gca;

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
% Larger font since this will be a sub panel
xl = xlabel('Temperature $\bar{T}$ [$^{\circ}$C]'); %textset(xl)
set(xl, 'Fontsize', 22, 'Interpreter', 'Latex');

yl = ylabel('$dV/dt$ [mm SLE yr$^-1$]'); %textset(yl)
set(yl, 'Fontsize', 22, 'Interpreter', 'Latex');
l1 = legend(num2str(legs, '%1.1f'), 'Location', 'SouthEast'); legset(l1)
legpos = l1.Position;
xlim([-1 4])

% Now do the picture-in-picture. Values are chosen for visibility.
xpos = 1;                 % Left dge
w = 2.75;                      % Width
ypos = -3.5;                % Lower edge
h = 4;                      % Height
% Draw rectangle around
rectangle('Position', [xpos ypos w h])%, 'LineStyle', '-.')
% Find T at middle of the box
midboxtemp = xpos + 0.5*w;

% Make lines from box to insert corners. First get corners of box
xboxupleft = xpos;
yboxupleft = ypos + h;
xboxdownright = xpos + w;
yboxdownright = ypos;

% Normalized figure units
[xbul, ybul] = ds2nfu(xboxupleft, yboxupleft);
[xbdr, ybdr] = ds2nfu(xboxdownright, yboxdownright);

% Make box for plot-in-plot; [left bottom width height]
%ax2 = axes('Position', [0.2, 0.15  0.45 0.45]);
insleft = 0.2;
insbottom = 0.15;
inswidth = 0.4;
insheight = 0.4;



% We are plotting this curve
interestvec = f(middle,:);
ax2 = axes('Position', [insleft insbottom  inswidth insheight]);
hold on; box on
plot(ax2, Tvec, interestvec, 'Color', keepcol, 'LineWidth', 2.5)
xlim([xpos, xpos+w])
ylim([ypos, ypos+h])

% Make lines from corner to corner
annotation('line', [xbul, insleft], [ybul, insbottom+insheight], ...
    'LineStyle', ':')
annotation('line', [xbdr, insleft+inswidth], [ybdr, insbottom], ...
    'LineStyle', ':')

%Plot middle, +/- sigma
[~, ind0] = min(abs(Tvec - midboxtemp));
[~, indP] = min(abs(Tvec - midboxtemp - sqrt(arvar)));
[~, indM] = min(abs(Tvec - midboxtemp + sqrt(arvar)));
% 
%hold on
scatter(Tvec(ind0), interestvec(ind0), 'ko')
scatter(Tvec(indP), interestvec(indP), 'ko')
scatter(Tvec(indM), interestvec(indM), 'ko')
% 
% Calculate f(x) = ax + b;
% a = (y2 - y1) / (x2 - x1)
ap = (interestvec(indP) - interestvec(indM)) / (Tvec(indP) - Tvec(indM));
bp = interestvec(indP) - ap * Tvec(indP);
fp = @(x) ap.*x + bp;
plot(Tvec, fp(Tvec), 'k-.')
% 
% %http://www.mathworks.com/matlabcentral/answers/4515-removing-ticks
set(gca, 'xtick', [])
set(gca, 'ytick', [])

% Arrows
% http://www.mathworks.com/matlabcentral/fileexchange/
% 10656-data-space-to-figure-units-conversion

yl = ylim;
ydist = yl(2) - yl(1);
sigstr = sprintf('$\\sigma_T$');%, 'Fontsize', fs, 'Interpreter', 'Latex');

xa1 = [midboxtemp - sqrt(arvar); midboxtemp-0.005];
yp1 = (yl(1) + 0.6*(yl(2) - yl(1)));
ya1 = interestvec(ind0)*[1 1];

[xaf1, yaf1] = ds2nfu(xa1, ya1);
arrow1 = annotation('doublearrow', xaf1, yaf1);
set(arrow1, 'LineWidth', 0.25)
text(midboxtemp - 0.5*sqrt(arvar), yp1-0.05*ydist, sigstr, 'Fontsize', fs, ...
   'Interpreter', 'Latex');

xa2 = [midboxtemp+0.03; midboxtemp + sqrt(arvar)];
yp2 = yl(1) + 0.5*(yl(2) - yl(1));
ya2 = interestvec(ind0)*[1 1];

[xaf2, yaf2] = ds2nfu(xa2, ya2);
arrow2 = annotation('doublearrow', xaf2, yaf2 );
text(midboxtemp + 0.5*sqrt(arvar), yp2+0.05*ydist, sigstr, 'Fontsize', fs, ...
   'Interpreter', 'Latex');
% 
% 
% Vertical dashed lines
l1 = line(Tvec(indM)*[1 1], [ya1(1) interestvec(indM)]);
set(l1, 'Color', 'k', 'LineStyle', '-.')
%
l2 = line(Tvec(indP)*[1 1], [ya1(2) interestvec(indP)]);
set(l2, 'Color', 'k', 'LineStyle', '-.')


% Write 'Volume [m. SLE]' over the legend box
%set('gca', orgplot)
text(orgplot, 2.7, -6.8, ...
    'Volume [m. SLE]', 'Fontsize', fs, 'Interpreter', 'Latex',...
    'BackgroundColor', 'w');

%% Save figures?
% figure names correspond to AGU figure naming convention
if save_pdf
    fprintf('Saving pdfs...\n');
    % pdf files for the article
    export_fig(fig010, [pdfpath, '2016gl070016-p02.pdf']);
    %print(fig011, [figpath, 'Sim+Approx-2016.eps'], '-depsc', '-r400')
    print(fig011, [pdfpath, '2016gl070016-p01.pdf'], '-dpdf', '-r400')
end

% png files used to generate README.md
if save_png
    fprintf('Saving pngs...\n');
    export_fig(fig010, [pngpath, '2016gl070016-p02.png']);
    % export_fig does not work with this figure
    print(fig011, [pngpath, '2016gl070016-p01.png'], '-dpng', '-r200')
    print(fig010, [pngpath, 'figure-two-left-panel.jpg'], '-djpeg', '-r600')
    export_fig(fig010, [pngpath, '2016gl070016-p02.png']);
end

display(temps)
display(meanV)

display(temps_show)
display(meansAR_show)

% -------------------------------------------------------------------------
fprintf('Done.\n')
toc
% -------------------------------------------------------------------------

