% figC1_wba_lyap.m
% WBA <-> Lyapunov cross-check for the standard map at K = 6.9.
% Two independent chaos indicators, computed on the same 120x120 = 14,400
% grid of initial conditions, are compared point by point:
%   x-axis: Benettin largest Lyapunov exponent lambda
%   y-axis: WBA convergence digits dig = -log10 |WBA[0,T]-WBA[T,2T]|
% Classification boundaries: dig >= 5 (regular, WBA); lambda < 0.03 (regular,
% Benettin). Points collapse into the two agreeing quadrants
% (both regular / both chaotic); agreement is 99.99%.
%
% Data: paper/data_matlab/sv3_wba_lyap.mat  (dig, lam, thresholds).
% This is a cross-check, not a discovery.

clear; close all;

here = fileparts(mfilename('fullpath'));
S = load(fullfile(here, '..', 'data_matlab', 'sv3_wba_lyap.mat'));

dig = double(S.dig(:));
lam = double(S.lam(:));            % 'lambda' is reserved in MATLAB
digcut = double(S.chaos_thresh);   % 5
lamcut = double(S.lyap_thresh);    % 0.03

wba_reg  = dig >= digcut;
lyap_reg = lam <  lamcut;
both_reg   = wba_reg & lyap_reg;
both_chaos = ~wba_reg & ~lyap_reg;
disagree   = wba_reg ~= lyap_reg;
agreement  = mean(wba_reg == lyap_reg);

% --- colours (colourblind-safe; no jet) ---
cReg   = [0.11 0.45 0.74];   % blue  - both regular
cChaos = [0.85 0.37 0.10];   % orange - both chaotic
cBad   = [0.15 0.15 0.15];   % dark  - disagreement

fig = figure('Units','centimeters','Position',[2 2 12 10],'Color','w');
ax = axes(fig); hold(ax,'on'); box(ax,'on');

% log-x helps because lambda spans ~1e-3 to ~1.4
set(ax,'XScale','log');

scatter(ax, max(lam(both_chaos),1e-4), dig(both_chaos), 6, cChaos, 'filled', ...
        'MarkerFaceAlpha',0.35,'MarkerEdgeColor','none');
scatter(ax, max(lam(both_reg),1e-4),   dig(both_reg),   10, cReg, 'filled', ...
        'MarkerFaceAlpha',0.85,'MarkerEdgeColor','none');
if any(disagree)
    scatter(ax, max(lam(disagree),1e-4), dig(disagree), 42, cBad, 'x', 'LineWidth',1.3);
end

% classification boundaries
yl = [-1 16]; xl = [8e-4 2];
plot(ax, xl, [digcut digcut], 'k--', 'LineWidth',0.8);
plot(ax, [lamcut lamcut], yl, 'k:',  'LineWidth',0.8);

xlim(ax, xl); ylim(ax, yl);
xlabel(ax, 'Benettin Lyapunov exponent  \lambda');
ylabel(ax, 'WBA convergence digits  dig_{1000}');
title(ax, 'WBA vs Lyapunov, standard map  K = 6.9');

% quadrant labels
text(ax, 1.2e-3, 13.5, 'both regular', 'Color',cReg,'FontWeight','bold','FontSize',9);
text(ax, 0.25, 1.2, 'both chaotic', 'Color',cChaos,'FontWeight','bold', ...
     'FontSize',9,'HorizontalAlignment','right');

nboth_reg = sum(both_reg); nboth_chaos = sum(both_chaos); nbad = sum(disagree);
msg = sprintf(['two independent chaos indicators agree\n' ...
               'agreement %.2f%% (14,400 ICs)\n' ...
               'both regular %d  /  both chaotic %d  /  disagree %d'], ...
               100*agreement, nboth_reg, nboth_chaos, nbad);
text(ax, 9e-4, -0.2, msg, 'FontSize',8,'VerticalAlignment','top', ...
     'BackgroundColor',[1 1 1 0.6]);

set(ax,'FontSize',9,'Layer','top','TickDir','out');

% many points -> raster the content at 300 dpi, keep vector axes/text
outpdf = fullfile(here, 'figC1_wba_lyap.pdf');
exportgraphics(fig, outpdf, 'ContentType','image', 'Resolution',300);
fprintf('wrote %s  (agreement %.4f)\n', outpdf, agreement);
