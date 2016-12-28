% 12 month average
% --> sigma
clearvars; close all; clc; format compact
rng('default')

% -------------------------------------------------------------------------
tic;
% -------------------------------------------------------------------------


% ------------- Save figures? ---------------------------------------------
savefigures = false; 
%
if savefigures == true;
    display('*--------------> Saving figures <--------------*')
end
% -------------------------------------------------------------------------

% ------------- Save data? ------------------------------------------------
savedata = true; 
%
if savedata == true;
    display('*--------------> Saving data <-----------------*')
end
% -------------------------------------------------------------------------

% ------------- Pretty figures --------------------------------------------
fs = 18;
figset = @(f)   set(f, 'Color', 'w');
textset = @(h)  set(h, 'Fontsize', fs, 'Interpreter', 'Latex');
legset = @(l)   set(l, 'FontSize', fs, 'Interpreter', 'Latex');
hcbset = @(h) set(h, 'FontSize', fs, 'FontName', 'Times New Roman');
set(0,'DefaultAxesFontSize',12)
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
addpath(genpath('MATLAB-Utilities'));
% Current working directory
currentwd = pwd;
% Parent folder
parentf = fileparts(currentwd);
% Where figures and data is stored
figpath = [parentf, '/LaTeX[IceSheets][2016]/gfx/'];
datapath = [parentf, '/Data[IceSheets][2016]/'];
figformat = '.pdf';
% Save results of ar1 fitting here
arma11filename = [datapath, 'arma11results.mat'];
% -------------------------------------------------------------------------


% ------------- Prettier colors -------------------------------------------
cls = lines(5);
blue   = cls(1, :);
red    = cls(2, :);
yellow = cls(3, :);
purple = cls(4, :);
green  = cls(5, :);
grey   = 0.4*[1 1 1]; 
% -------------------------------------------------------------------------






%% Scan file
fileid = fopen([datapath, 'GreenlandTemperatureData.txt']);
C = textscan(fileid, '%f %f', 'CommentStyle', '#');
fclose(fileid);

% Year and temp
year = C{1};
fyear = floor(year);
temp = C{2};


% Want the average temerature each year;
% 12 month average 
uyear = unique(fyear);
avgtemp = zeros(size(uyear));
for ii = 1:numel(uyear)
    temps = temp(fyear == uyear(ii));
    avgtemp(ii) = mean(temps);
end

% Standard deviation
stdtemp = std(avgtemp);

% Figure: monthly temp and yearly means
fig001 = figure(001); hold on; box on; figset(fig001)
p1 = plot(year, temp, 'color', blue);
p2 = plot(uyear, avgtemp, 'color', red, 'linewidth', 2);
% Pretty
%tt = title('Observed temperature anomaly over Greenland'); textset(tt)
xl = xlabel('Year'); textset(xl)
yl = ylabel('Temperature Anomaly [$^{\circ}$C]'); textset(yl)
l1 = legend([p1 p2], {'Monthly mean', 'Yearly mean'}); legset(l1)
xlim([uyear(1) uyear(end)])



%% Estimate autocorrelation of Greenland temperature
nlags = floor(numel(avgtemp)/4);
[acf, lags, bounds] = autocorr(avgtemp, nlags);


% Figure: ACF of temperature
fig3 = figure(3); hold on; box on; figset(fig3)
sc = scatter(lags, acf, 'ko');
oe = repmat(exp(-1), nlags+1, 1);
%poe = plot(lags, oe);


% Find the first temp below 1/e
cct = find(acf < exp(-1), 1) - 1; % Minus 1 due to lag 0 being included
lexp = line(xlim, [exp(-1) exp(-1)], 'Color', red);
lcct = line([cct cct], ylim, 'Color', blue);

legstrs = {'ACF$(T_t)$', '$1/e$', 'CCT$(T_t)$'};
l1 = legend([sc lexp lcct], legstrs); legset(l1)

tt = title('Autocorrelation of observed temperatures'); textset(tt)
xl = xlabel('Lag $k$'); textset(xl)
yl = ylabel('ACF$(k)$'); textset(yl)

%% Fit ARMA(1, 1) process
armamdl = arima(1, 0, 1);
% Estimate model. len(v) must be =(<=?) len(t) - 1
armaest = estimate(armamdl, avgtemp, 'Display', 'off');
% Get the parameters. Note the beautiful syntax of MATLAB
estar1 = armaest.AR{1};
estma1 = armaest.MA{1};
estvar = armaest.Variance;

% Save the varince and ar1-parameter
arma11 = struct;
arma11.estar1 = estar1;
arma11.estma1 = estma1;
arma11.estvar = estvar;

% Save the results
if savedata == true
    % ARMA parameters
    save(arma11filename, '-struct', 'arma11');
end


%% Simulate artificial temperatures
tsim = zeros(numel(avgtemp), 1);
tsim(1) = 0;
sqrtvar = sqrt(estvar);
eta = sqrtvar*randn(numel(avgtemp)+1);


% Note that the noise vector has been extended to length numel(avgtemp)+1
% so that in the noise vector "now" is eta(idx+1) and "one step back" is
% eta(idx)


% Manual AR1
for idx = 2:numel(avgtemp)
    tsim(idx) = ...
        estar1*tsim(idx-1) + ...    % AR part of model
        estma1*eta(idx) + ...       % MA part of model
        eta(idx+1);                 % White nise
end

% Figure: observed and simulated temperature
fig004 = figure(004); hold on; box on; figset(fig004)
pobs = plot(uyear, avgtemp, 'color', grey, 'linewidth', 2);
psim1 = plot(uyear, tsim(:, 1), 'color', blue);
% Pretty plot
legstrs = {'Observed', 'Simulated ARMA(1, 1) process'};
l1 = legend([pobs psim1], legstrs); legset(l1)
%tt = title('Simulation of Greenland temperature anomaly'); textset(tt)
xl = xlabel('Year'); textset(xl)
yl = ylabel('Temperature Anomaly [$^{\circ}$C]'); textset(yl)

%% Variance 
%  http://faculty.chicagobooth.edu/ruey.tsay/teaching/uts/lec2-08.pdf, p. 7
%  note the sign change in the numerator, compare with definition of theta
%  https://se.mathworks.com/help/econ/arima-class.html under
%  "Linear Time Series Model"
armavar = @(sig, phi, theta) sig^2 * ...
    (1 + theta^2 + 2*phi*theta) / (1 - phi^2);

% Compute the theoretical variance of the ARMA(1, 1) process
theovar = armavar(sqrtvar, estar1, estma1);

% Generate nsim temperature series, compute the variance, compare
nsim = 10;
simlen = 5000;
simvars = zeros(nsim, 1);

for idy = 1:nsim
    
    tsim = zeros(simlen, 1);
    tsim(1) = 0;
    eta = sqrtvar*randn(simlen + 1);

    % Manual ARMA(1, 1)
    for idx = 2:simlen
        tsim(idx) = ...
            estar1*tsim(idx-1) + ...    % AR part of model
            estma1*eta(idx) + ...       % MA part of model
            eta(idx+1);                 % White nise
    end
    
    simvars(idy) = var(tsim);
    
end


%%
fig005 = figure(005); figset(fig005); hold on; box on

sc = scatter(1:1:nsim, simvars, 'ko');
ylim([0 2])

ltheo = line(xlim, [theovar theovar], 'color', red);
lsim = line(xlim, [mean(simvars) mean(simvars)], 'color', green);

legstrs = {'$\sigma^2 \left( T_t^{(i)} \right)$', ...
    'Theoretical $\sigma^2 \left( T_t \right)$', ...
    '$\mu \left[ \sigma^2 \left( T_t^{(i)} \right) \right]$'};
l1 = legend([sc ltheo lsim], legstrs, 'location', 'southeast'); legset(l1)

xl = xlabel('Simulation index $i$'); textset(xl)
yl = ylabel('Variance'); textset(yl)


%% ACF(Tt)
simlen = 5e2;
tsim = zeros(simlen, 1);
tsim(1) = 0;
sqrtvar = sqrt(estvar);
eta = sqrtvar*randn(numel(avgtemp)+1);

% Manual ARMA(1, 1)
for idx = 2:simlen
    tsim(idx) = ...
        estar1*tsim(idx-1) + ...    % AR part of model
        estma1*eta(idx) + ...       % MA part of model
        eta(idx+1);                 % White nise
end



% Estimate autocorrelation of Greenland temperature
nlags = floor(numel(avgtemp)/4);
[acf, lags, bounds] = autocorr(tsim, nlags);


% Figure: ACF of temperature
fig006 = figure(006); hold on; box on; figset(fig006)
sc = scatter(lags, acf, 'ko');
oe = repmat(exp(-1), nlags+1, 1);
%poe = plot(lags, oe);


% Find the first temp below 1/e
cct = find(acf < exp(-1), 1) - 1; % Minus 1 due to lag 0 being included
lexp = line(xlim, [exp(-1) exp(-1)], 'Color', red);
lcct = line([cct cct], ylim, 'Color', blue);

legstrs = {'ACF$(T_t)$', '$1/e$', 'CCT$(T_t)$'};
l1 = legend([sc lexp lcct], legstrs); legset(l1)

tt = title('Autocorrelation of simulated ARMA(1, 1) temperatures'); textset(tt)
xl = xlabel('Lag $k$'); textset(xl)
yl = ylabel('ACF$(k)$'); textset(yl)
