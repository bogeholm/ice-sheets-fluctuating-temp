% 12 month average
% --> sigma
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
% Save results of ar1 fitting here
ar1filename = [datapath, 'ar1results.mat'];
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
savefigures = true;
%savefigures = false;
%savedata = true;
savedata = true;
% 'test' is quite quick (~20 seconds), 'production' around 5 minutes
thisrun = 'test';
%thisrun = 'production';
% -------------------------------------------------------------------------



% ------------- Verbosity ---------------------------------------------
if savefigures == true
    disp('*--------------> Saving figures <--------------*')
end


if savedata == true
    disp('*--------------> Saving data <-----------------*')
end
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


% Save temperatures for use in R
%if savedata == true
%    csvwrite([datapath, 'greenlandAVGtemp.csv'], [uyear avgtemp])
%end


%% Estimate autocorrelation of Greenland temperature
nlags = floor(numel(avgtemp)/4);
[acf, lags, bounds] = autocorr(avgtemp, nlags);


% Figure: ACF of temperature
fig3 = figure(3); hold on; box on; figset(fig3)
plot(lags, acf)
oe = repmat(exp(-1), nlags+1, 1);
%poe = plot(lags, oe);


% Find the first temp below 1/e
cct = find(acf < exp(-1), 1);
line([cct cct], ylim);


%% Fit AR(1) process to temperature
% Specify AR part of the model
mdl = arima(1, 0, 0);
% Estimate model. len(v) must be =(<=?) len(t) - 1
est = estimate(mdl, avgtemp, 'Display', 'off');
clc
% Get the parameters. Note the beautiful syntax of MATLAB
estar1 = est.AR{1};
estvar = est.Variance;

% Save the varince and ar1-parameter
ar1 = struct;
ar1.estar1 = estar1;
ar1.estvar = estvar;

% Save the results
if savedata == true
    % AR parameters
    save(ar1filename, '-struct', 'ar1');
end

%% Fit ARMA(1, 1) process
armamdl = arima(1, 0, 1);
% Estimate model. len(v) must be =(<=?) len(t) - 1
armaest = estimate(armamdl, avgtemp, 'Display', 'off');
% Get the parameters. Note the beautiful syntax of MATLAB
%estar1 = est.AR{1};
%estvar = est.Variance;

% Save the varince and ar1-parameter
%ar1 = struct;
%ar1.estar1 = estar1;
%ar1.estvar = estvar;

% Save the results
%if savedata == true
%    % AR parameters
%    save(ar1filename, '-struct', 'ar1');
%end


%%

% Simulate artificial temperatures
tsim = zeros(numel(avgtemp), 1);
tsim(1) = 0;
eta = randn(size(avgtemp));

% Manual AR1
for idx = 2:numel(avgtemp)
    tsim(idx) = estar1*tsim(idx-1) + sqrt(estvar)*eta(idx);
end

% Figure: observed and simulated temperature
fig004 = figure(004); hold on; box on; figset(fig004)
pobs = plot(uyear, avgtemp, 'color', grey, 'linewidth', 2);
psim1 = plot(uyear, tsim(:, 1), 'color', blue);
% Pretty plot
legstrs = {'Observed', 'Simulated AR(1) process'};
l1 = legend([pobs psim1], legstrs); legset(l1)
%tt = title('Simulation of Greenland temperature anomaly'); textset(tt)
xl = xlabel('Year'); textset(xl)
yl = ylabel('Temperature Anomaly [$^{\circ}$C]'); textset(yl)





%% Relationship between estvar and variance of the generated temperature

% ------------- Variance of AR(1) process ---------------------------------
%https://en.wikipedia.org/wiki/Autoregressive_model#Example:_An_AR.281.29_process
%
% x(t+1) = a1*x(t) + s*eta => var(xt) = s^2/(1 - a1^2)
arvar_func = @(a1, s) s^2 / (1 - a1^2);
% -------------------------------------------------------------------------
arvar = arvar_func(estar1, sqrt(estvar));


% Simulate artificial temperatures
nsims = 20;
lensim = 2500;
burnin = 500; % Let the process reach steady state
numericalvar = zeros(nsims, 1);

for simid = 1:nsims
    tsim_check=zeros(lensim, 1);
    eta = randn(lensim);
    for idx = 2:lensim
        tsim_check(idx) = estar1*tsim_check(idx-1) + sqrt(estvar)*eta(idx);
    end
    numericalvar(simid) = var(tsim_check(burnin:end));
end

mvar = mean(numericalvar);


%{
% Plot variance of simulated temperatures
xax = 1:1:nsims;
one = ones(size(xax));

% Figure: comparison of variances
fig5 = figure(5); hold on; box on; figset(fig5)

pl1 = plot(xax, one*arvar, 'color', blue, 'linewidth', 1.5);
pl2 = plot(xax, one*mvar, 'color', green, 'linewidth', 1.5);
pl3 = scatter(xax, numericalvar, 'markerfacecolor', red);

% Pretty plot
legstrs = {'Theory', 'Mean Variance of Simulations', ...
    'Variance of Simulations'};
l1 = legend([pl1, pl2, pl3], legstrs); legset(l1)
tt = title('Variance of Temperature Simulations'); textset(tt)
xl = xlabel('Simulation Number'); textset(xl)
yl = ylabel('Variance of $T$'); textset(yl)
ylim([1 2])
%}

%% Save figures?
if savefigures
    fprintf('Saving figures...\n');
    %export_fig(fig001, [figpath, 'GreenlandTemp-2016.pdf'])
    export_fig(fig001, [figpath, 'GreenlandTemp-2016.png'])
end





%% Print out data
fprintf('\n\n    ****************** Computed Values *****************\n\n')

fprintf('----------- Variance of Yearly Mean Anomalies ----- \n')
fprintf('Var(T): %.4f', var(avgtemp));


fprintf('\n\n----------- AR(1)-model ---------------------------\n')
fprintf('Model is: \n');
fprintf('\t 1) T_{t+1} = phi*T_t + S*W_t\n');
fprintf('\n');
fprintf('Estimated parameters: \n');
fprintf('phi = %.4f\n', estar1);
fprintf('S^2 = %.4f\n', estvar);
fprintf('\n')
fprintf('Theoretical variance of AR(1) model:\n')
fprintf('Var(theory) = %.4f\n', arvar);

fprintf('\n\n----------- Mean variance of %i {T_t} --------------\n', nsims)
fprintf('            series simulated with 1)\n')
fprintf('Var({T_t}): %.4f\n', mean(numericalvar));


% -------------------------------------------------------------------------
fprintf('Done.\n');
toc
% -------------------------------------------------------------------------

