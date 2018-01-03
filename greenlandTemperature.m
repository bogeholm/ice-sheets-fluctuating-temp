%% Load Greenland temperature data and fit AR(1)-model
%  to annual averages
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


% -------------------------------------------------------------------------
% Save results of ar1 fitting here
ar1filename = [datapath, 'ar1results.mat'];
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
savedata = true;
%savedata = false;
% -------------------------------------------------------------------------
hlred_rgb = [174 0 20];
hlred = hlred_rgb / 255;


% ------------- Verbosity ---------------------------------------------
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
p2 = plot(uyear, avgtemp, 'color', hlred, 'linewidth', 2);
% Pretty
%tt = title('Observed temperature anomaly over Greenland'); textset(tt)
xl = xlabel('Year'); textset(xl)
yl = ylabel('Temperature Anomaly ($^{\circ}$C)'); textset(yl)
l1 = legend([p1 p2], {'Monthly mean', 'Annual mean'}); legset(l1)
xlim([uyear(1) uyear(end)])




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
fig002 = figure(002); hold on; box on; figset(fig002)
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




% Save figures?
if save_pdf
    fprintf('Saving pdfs...\n');
    export_fig(fig001, [pdfpath, 'figureS02.pdf'])
end

% Save figures?
if save_png
    fprintf('Saving pngs...\n');
    export_fig(fig001, [pngpath, 'figureS02.png'])
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




% -------------------------------------------------------------------------
fprintf('Done.\n');
toc
% -------------------------------------------------------------------------

