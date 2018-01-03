%% Check the consequences of the time step size used to 
%  integrate the Oerlemans model
%
% Troels B. Mikkelsen - bogeholm@nbi.ku.dk
% September 2016


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



% Mean temperature for integration
temps = 0;
% Step sizes 
stepsizes = [1; 0.5; 0.1; 0.01];
% End of integration
timemax      = 1e4;
n_temps     = numel(temps); % This many temperatures
% -------------------------------------------------------------------------



% ------------- Variance of AR(1) process ---------------------------------
%https://en.wikipedia.org/wiki/Autoregressive_model#
% Example:_An_AR.281.29_process
%
% x(t+1) = a1*x(t) + s*eta => var(xt) = s^2/(1 - a1^2)
arvar_func = @(a1, s) s^2 / (1 - a1^2);
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Volume as function of radius; 
Vr = @(R) iceTotal(R, par);
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


Rstart = 7e5;   % Initial volume; inspect results to ensure sanity
Tsim = 0;

VolAndTime = cell(numel(stepsizes), 2);

% Run the forward Euler on O/03
for idstep = 1:numel(stepsizes)
    dt = stepsizes(idstep);
    Rval    = Rstart; % initial conditions
    time    = -dt;
    ii      = 0;
    %
    timeax = (0:dt:timemax)';
    steps = numel(timeax);
    
    %
    Rresults = nan(numel(timeax), 1);
    % Do forward euler
    while time < timemax
        time    = time + dt;
        ii      = ii + 1;
        [dRdt, dVdt] = oerlemansModel(time, Rval, Tsim, par);
        Rval = Rval + dt*dRdt;
        % Store
        Rresults(ii) = Rval;
     end % end while
     % O/03 works in radius - we convert
    Vresults = OceanSurf*Vr(Rresults);
    % Save this
    VolAndTime{idstep, 1} = Vresults; 
    VolAndTime{idstep, 2} = timeax; 
end


%% Plot this
fig001 = figure(001); hold on; box on; figset(fig001)

for idstep = 1:numel(stepsizes)
    v = VolAndTime{idstep, 1};
    t = VolAndTime{idstep, 2};
    plot(t, v)
end


% Legend strings
legstr = cell(numel(stepsizes),1);
for idx = 1:numel(stepsizes)
    str = ['$\Delta t = ', sprintf('%.2f', stepsizes(idx)), '$'];
    legstr{idx} = str;
end
lh = legend(legstr, 'location', 'northwest'); legset(lh)
xl = xlabel('Time [yr]'); textset(xl)
yl = ylabel('Volume [m. SLE]'); textset(yl)



%% Vary the temperature ever so gently
disp('Integrating slowly varying temperatures')
% We use a sine function
P = 500; % Period
A = 1; % Amplitude
Tfunc = @(t) temps + A*sin(2*pi*t/P);


VolAndTime_sin = cell(numel(stepsizes), 2);

% Run the forward Euler on O/03
for idstep = 1:numel(stepsizes)
    dt = stepsizes(idstep);
    Rval    = Rstart; % initial conditions
    time    = -dt;
    ii      = 0;
    
    % disp(dt)
    
    timeax = (0:dt:timemax)';
    steps = numel(timeax);
    Tsim = Tfunc(timeax);
    %
    Rresults = nan(numel(timeax), 1);
    % Do forward euler
    while time < timemax
        time    = time + dt;
        ii      = ii + 1;
        [dRdt, dVdt] = oerlemansModel(time, Rval, Tsim(ii), par);
        Rval = Rval + dt*dRdt;
        % Store
        Rresults(ii) = Rval;
     end % end while
     % O/03 works in radius - we convert
    Vresults = OceanSurf*Vr(Rresults);
    % Save this
    VolAndTime_sin{idstep, 1} = Vresults; 
    VolAndTime_sin{idstep, 2} = timeax; 
end



%% Plot slowly varying temperatures
fig002 = figure(002); hold on; box on; figset(fig002)

for idstep = 1:numel(stepsizes)
    v = VolAndTime_sin{idstep, 1};
    t = VolAndTime_sin{idstep, 2};
    plot(t, v)
end


% Legend strings
legstr = cell(numel(stepsizes),1);
for idx = 1:numel(stepsizes)
    str = ['$\Delta t = ', sprintf('%.2f', stepsizes(idx)), '$'];
    legstr{idx} = str;
end
lh = legend(legstr, 'location', 'northeast'); legset(lh)
xl = xlabel('Time [yr]'); textset(xl)
yl = ylabel('Volume [m. SLE]'); textset(yl)

%% Now, generate AR1-noise and use these results
%  1) Generate AR(1) time series, one temp per year
%  2) Pad this so it matches the time step size - ie. 100 consecutive temps
%  for dt = 0.01
%  3) Integrate!
disp('Generating AR(1) temperatures')

meantemp = temps;

% Simulate AR(1) temperature
% Generate for one extra year
tempar = zeros(timemax+1, 1);
eta = randn(timemax+1, 1); % noise
c = meantemp*(1 - estar1);
tempar(1) = meantemp;
% AR(1) temperature tie series
for idx = 2:numel(tempar)
    tempar(idx) = c + estar1*tempar(idx-1) + sqrt(estvar)*eta(idx);
end


% Store results here
TempsAR = cell(numel(stepsizes), 1);


for idstep = 1:numel(stepsizes)
    dt = stepsizes(idstep);
    N = 1/dt;
    % disp(N)
    % https://se.mathworks.com/matlabcentral/answers/
    % 46898-repeat-element-of-a-vector-n-times-without-loop
    there = repmat(tempar, 1, N)';
    there = there(:)';

    TempsAR{idstep, 1} = there;
end



% Plot temperatures, the start only
% fig003 = figure(003); hold on; box on; figset(fig003)
% 
% plotlim = 10; % plot until plotlim years
% 
% for idstep = 1:numel(stepsizes)
%     dt = stepsizes(idstep);
%     
%     t = 0:dt:plotlim;
%     
%     there = TempsAR{idstep, 1};
%     plot(t, there(1:numel(t)))
% end



%% Use these temps to integrate
disp('Integrating AR(1) temperatures')
VolAndTime_AR = cell(numel(stepsizes), 2);

% Run the forward Euler on O/03
for idstep = 1:numel(stepsizes)
    dt = stepsizes(idstep);
    Rval    = Rstart; % initial conditions
    time    = -dt;
    ii      = 0;
    
    %disp(dt)
    
    timeax = (0:dt:timemax)';
    %steps = numel(timeax);
    there = TempsAR{idstep, 1};
    %
    Rresults = nan(numel(timeax), 1);
    % Do forward euler
    while time < timemax
        time    = time + dt;
        ii      = ii + 1;
        [dRdt, dVdt] = oerlemansModel(time, Rval, there(ii), par);
        Rval = Rval + dt*dRdt;
        % Store
        Rresults(ii) = Rval;
     end % end while
     % O/03 works in radius - we convert
    Vresults = OceanSurf*Vr(Rresults);
    % Save this
    VolAndTime_AR{idstep, 1} = Vresults; 
    VolAndTime_AR{idstep, 2} = timeax; 
end



%% Plot this
fig004 = figure(004); hold on; box on; figset(fig004)

for idstep = 1:numel(stepsizes)
    v = VolAndTime_AR{idstep, 1};
    t = VolAndTime_AR{idstep, 2};
    
    plot(t, v)
end


% Legend strings
legstr = cell(numel(stepsizes),1);
for idx = 1:numel(stepsizes)
    str = ['$\Delta t = ', sprintf('%.2f', stepsizes(idx)), '$'];
    legstr{idx} = str;
end

lh = legend(legstr, 'location', 'northeast'); legset(lh)
xl = xlabel('Time [yr]'); textset(xl)
yl = ylabel('Volume [m. SLE]'); textset(yl)


%% Calculate the mean absolute deviation from dt = 1 year
disp('Calculating mean absolute deviation')
v1 = VolAndTime_AR{1, 1};
meanv = mean(v1);


for idstep = 2:numel(stepsizes)
    dt = stepsizes(idstep);
    N = 1/dt;
    
    % Get volume
    vhere = VolAndTime_AR{idstep, 1};
    
    % repmat v1 to match number of elements in vhere
    vone = repmat(v1, 1, N)';
    vone = vone(:)';
    % Volume for the last time step is repeated one too many times
    vone = vone(1:end-(N-1));
    
    % Column or row?
    if ~all(size(vone) == size(vhere))
        vone = vone';
    end
    
    % Absolute difference
    absdiff = abs(vhere - vone);
    
    % Verbosity
    fprintf('\n')
    fprintf('dt = %.2f\n', dt)
    disp('mean(abs(diff))')
    disp(mean(absdiff))
    disp('mean(abs(diff))/mean(V(dt=1))')
    disp(mean(absdiff)/meanv)
    fprintf('\n')
end






%% Zoom in on a random part of the above plot
span = 15; % Number of years we will consider
plotyear = randi(timemax - span);

fig005 = figure(005); figset(fig005)
subplot(211); hold on; box on
subplot(212); hold on; box on

% Markersizes. 6 is default.
marks = [6, 5, 4, 1];

for idstep = 1:numel(stepsizes)
    dt = stepsizes(idstep);
    
    vv = VolAndTime_AR{idstep, 1};
    tt = VolAndTime_AR{idstep, 2};
    T = TempsAR{idstep, 1};
    
    % We only plot these
    mrbool = tt >= plotyear & tt <= plotyear+span;
    
    % Volume
    subplot(211)
    plot(tt(mrbool), vv(mrbool), '-o', 'MarkerSize', marks(idstep))
    
    % Temperature
    subplot(212)
    plot(tt(mrbool), T(mrbool))
end



% Legend strings
legstr = cell(numel(stepsizes),1);
for idx = 1:numel(stepsizes)
    str = ['$\Delta t = ', sprintf('%.2f', stepsizes(idx)), '$'];
    legstr{idx} = str;
end

subplot(211)
lh = legend(legstr, 'location', 'northwest'); legset(lh)
xl = xlabel('Time [yr]'); textset(xl)
yl = ylabel('Volume [m. SLE]'); textset(yl)

subplot(212)
lh = legend(legstr, 'location', 'northwest'); legset(lh)
xl = xlabel('Time [yr]'); textset(xl)
yl = ylabel('Temp [$^{\circ}$C]'); textset(yl)





%% One more zoom - now only dt=1 and dt=0.01
fig006 = figure(006); figset(fig006); hold on; box on

% Plot from start to stop
start = 0;
stop = 200;

% Time/Volume for dt = 1
v1 = VolAndTime_AR{1, 1};
t1 = VolAndTime_AR{1, 2};
% Time/Volume for dt = 0.001
v001 = VolAndTime_AR{4, 1};
t001 = VolAndTime_AR{4, 2};
% Relevant indices
mrbool1 = t1 >= start & t1 <= stop;
mrbool001 = t001 >= start & t001 <= stop;

% Plot
plot(t001(mrbool001), v001(mrbool001), 'Color', red, 'LineWidth', 2)
sc = scatter(t1(mrbool1), v1(mrbool1), 'MarkerEdgeColor', blue);


lh = legend('$\Delta t = 1$', '$\Delta t = 0.01$', ...
    'location', 'northeast'); legset(lh)

xl = xlabel('Time [yr]'); textset(xl)
yl = ylabel('Volume [m. SLE]'); textset(yl)



% ------------------------------------------------------------------------
if save_pdf
    disp('Saving pdfs...')
    export_fig(fig004, [pdfpath, 'figureS01.pdf']);
end
% ------------------------------------------------------------------------
if save_png
    disp('Saving pngs...')
    export_fig(fig004, [pngpath, 'figureS01.png']);
end
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
toc;
% -------------------------------------------------------------------------




