% ------------- Pretty figures --------------------------------------------
fs = 14;
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
grey = 0.4*[1, 1, 1];
lightgrey   = 0.25*[1, 1, 1];
darkgrey   = 0.55*[1, 1, 1];
% http://cloford.com/resources/colours/500col.htm
violetred = [255, 62, 150] / 255; % violetred 1
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Where figures and data is stored
figpath = 'gfx/';
texpath = 'tex/';
datapath = 'data/';
robinsonpath = 'robinsondata/';
figformat = '.pdf';
% -------------------------------------------------------------------------