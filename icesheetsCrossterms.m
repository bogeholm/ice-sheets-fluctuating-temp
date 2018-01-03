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
%thisrun = 'test';
thisrun = 'production';
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
    disp('No file ar1results.mat - run greenlandTemperature.m')
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
yl = ylabel('$\left(V(t) - \overline{V}\right)^2$'); textset(yl)
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
yl = ylabel('$\left(V(t) - \overline{V}\right) \left(T(t) - \overline{T}\right)$'); textset(yl)
tt = title('$(V - \overline{V})(T - \overline{T})$ as a function of time'); 
textset(tt)
ylim([-20 20])
xlim([0 2e4])





%% Save pdf?
if save_pdf
    disp('Saving pdfs...')
    export_fig(fig117, [pdfpath, 'figureS03.pdf'])
end

if save_png
    disp('Saving pngs...')
    print(fig117, [pngpath, 'figureS03.png'], '-dpng', '-r600')

end

% -------------------------------------------------------------------------
fprintf('Done.\n')
toc
% -------------------------------------------------------------------------