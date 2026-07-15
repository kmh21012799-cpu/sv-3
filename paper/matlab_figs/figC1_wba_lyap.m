% figC1_wba_lyap.m
% WBA <-> Lyapunov cross-check for the standard map at K = 6.9.
% 14,400 initial conditions (120x120 grid); each point compared by two
% independent chaos indicators:
%   x-axis: Benettin largest Lyapunov exponent lambda
%   y-axis: WBA convergence digits dig = -log10 |WBA[0,T]-WBA[T,2T]|
% Classification boundaries: dig = 5 (WBA cut); lambda = 0.03 (Benettin cut).
% Points collapse into the two agreeing quadrants; agreement is 99.99%.
% Counts and the summary sentence go in the figure CAPTION, not on the axes.
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
cBad   = [0.10 0.10 0.10];   % black  - disagreement
cGrid  = [0.45 0.45 0.45];   % grey   - boundary-line labels

fig = figure('Units','centimeters','Position',[2 2 12 9.5],'Color','w');
ax = axes(fig); hold(ax,'on'); box(ax,'on');
set(ax,'XScale','log');       % lambda spans ~1e-3 to ~1.4

% dense chaotic cloud: small marker + low alpha so it does not smear
scatter(ax, max(lam(both_chaos),1e-4), dig(both_chaos), 5, cChaos, 'filled', ...
        'MarkerFaceAlpha',0.22,'MarkerEdgeColor','none');
scatter(ax, max(lam(both_reg),1e-4),   dig(both_reg),   16, cReg, 'filled', ...
        'MarkerFaceAlpha',0.90,'MarkerEdgeColor','none');
if any(disagree)
    scatter(ax, max(lam(disagree),1e-4), dig(disagree), 46, cBad, 'x', 'LineWidth',1.4);
end

% classification boundaries (labelled)
yl = [-1 16]; xl = [8e-4 2];
plot(ax, xl, [digcut digcut], '--', 'Color',cGrid, 'LineWidth',0.7, 'HandleVisibility','off');
plot(ax, [lamcut lamcut], yl, ':',  'Color',cGrid, 'LineWidth',0.9, 'HandleVisibility','off');

xlim(ax, xl); ylim(ax, yl);
xlabel(ax, '\lambda   (Benettin Lyapunov exponent)');
ylabel(ax, 'dig_{1000}   (WBA convergence digits)');
set(ax,'FontSize',10,'Layer','top','TickDir','out');
set(ax,'XTick',[1e-3 1e-2 1e-1 1e0]);

% --- in-axes text: cluster labels + one agreement line + boundary labels only ---
text(ax, 1.15e-3, 14.6, 'both regular',  'Color',cReg,  'FontWeight','bold','FontSize',10);
text(ax, 1.75,    2.4,  'both chaotic',  'Color',cChaos,'FontWeight','bold','FontSize',10, ...
     'HorizontalAlignment','right');
text(ax, 1.75,   14.6,  '99.99% agreement', 'Color',[0.2 0.2 0.2],'FontSize',10, ...
     'HorizontalAlignment','right');
text(ax, 1.75,    6.4,  'disagree', 'Color',cBad,'FontSize',9, ...
     'HorizontalAlignment','right');
% boundary-line labels
text(ax, 9.2e-4, 5.5, 'dig = 5', 'Color',cGrid,'FontSize',8.5,'FontAngle','italic', ...
     'VerticalAlignment','bottom');
text(ax, 0.034, 8.5, '\lambda = 0.03', 'Color',cGrid,'FontSize',8.5,'FontAngle','italic');

% many points -> raster the content at 300 dpi, keep vector axes/text
outpdf = fullfile(here, 'figC1_wba_lyap.pdf');
exportgraphics(fig, outpdf, 'ContentType','image', 'Resolution',300);
fprintf('wrote %s  (agreement %.4f)\n', outpdf, agreement);
