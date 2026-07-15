% figC1_wba_lyap.m
% WBA <-> Lyapunov cross-check for the standard map at K = 6.9.
% Two independent chaos indicators, computed on the same 120x120 = 14,400
% grid of initial conditions, are compared point by point:
%   x-axis: Benettin largest Lyapunov exponent lambda
%   y-axis: WBA convergence digits dig = -log10 |WBA[0,T]-WBA[T,2T]|
% Classification boundaries: dig >= 5 (regular, WBA); lambda < 0.03 (regular,
% Benettin). Points collapse into the two agreeing quadrants
% (both regular / both chaotic); agreement is 99.99%.
% All annotation is outside the axes (legend at right).
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
cReg   = [0.11 0.45 0.74];   % blue   - both regular
cChaos = [0.85 0.37 0.10];   % orange - both chaotic
cBad   = [0.15 0.15 0.15];   % dark   - disagreement

fig = figure('Units','centimeters','Position',[2 2 14 9],'Color','w');
ax = axes(fig); hold(ax,'on'); box(ax,'on');
set(ax,'XScale','log');       % lambda spans ~1e-3 to ~1.4

h1 = scatter(ax, max(lam(both_chaos),1e-4), dig(both_chaos), 6, cChaos, 'filled', ...
        'MarkerFaceAlpha',0.35,'MarkerEdgeColor','none', ...
        'DisplayName',sprintf('both chaotic (%d)', sum(both_chaos)));
h2 = scatter(ax, max(lam(both_reg),1e-4),   dig(both_reg),   14, cReg, 'filled', ...
        'MarkerFaceAlpha',0.90,'MarkerEdgeColor','none', ...
        'DisplayName',sprintf('both regular (%d)', sum(both_reg)));
hs = [h1 h2];
if any(disagree)
    h3 = scatter(ax, max(lam(disagree),1e-4), dig(disagree), 42, cBad, 'x', ...
        'LineWidth',1.3, 'DisplayName',sprintf('disagree (%d)', sum(disagree)));
    hs = [h1 h2 h3];
end

% classification boundaries (kept out of the legend)
yl = [-1 16]; xl = [8e-4 2];
plot(ax, xl, [digcut digcut], 'k--', 'LineWidth',0.7, 'HandleVisibility','off');
plot(ax, [lamcut lamcut], yl, 'k:',  'LineWidth',0.7, 'HandleVisibility','off');

xlim(ax, xl); ylim(ax, yl);
xlabel(ax, '\lambda   (Benettin Lyapunov exponent)');
ylabel(ax, 'dig_{1000}   (WBA convergence digits)');
set(ax,'FontSize',10,'Layer','top','TickDir','out');
set(ax,'XTick',[1e-3 1e-2 1e-1 1e0]);

% all text outside the axes: legend at right, agreement as its title
lgd = legend(hs, 'Location','eastoutside');
lgd.Title.String = sprintf('agreement %.2f%%', 100*agreement);
lgd.FontSize = 9;

% many points -> raster the content at 300 dpi, keep vector axes/text
outpdf = fullfile(here, 'figC1_wba_lyap.pdf');
exportgraphics(fig, outpdf, 'ContentType','image', 'Resolution',300);
fprintf('wrote %s  (agreement %.4f)\n', outpdf, agreement);
