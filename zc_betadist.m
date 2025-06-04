%% * Preamble
clc
clear all
addpath('./ROUTINES/export_fig/')
addpath('./ROUTINES/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

savfigs = false;

%% * Parameter List
% alfs = [.5 3.5 2 5 .2 .5 .2 1];
% bts = [.5 3.5 5 2 .2 .2 .5 1];
alfs = [1 3.5 .5 5];
bts = [1 3.5 .5 2];

Nx = 100;
xs = linspace(0, 1, Nx)';
bps = zeros(Nx, length(alfs));
for bi=1:length(alfs)
    bps(:,bi) = betapdf(xs, alfs(bi), bts(bi));
end

%% * Plot
colos = DISTINGUISHABLE_COLORS(length(alfs))*0.85;

fsz = 18;
figure(1)
poss=get(gcf, 'Position');
set(gcf, 'Color', 'white', 'Position', poss)
clf()
tiledlayout(1,1,'TileSpacing','compact', 'Padding','compact');
nexttile;

aa = gobjects(length(alfs),1);
for bi=1:length(alfs)
    aa(bi) = plot(2*xs-1, bps(:, bi)/2, 'LineWidth', 2, 'Color', colos(bi,:)); hold on
    legend(aa(bi), sprintf('C%d: $(\\alpha, \\beta)=(%.1f, %.1f)$', ...
                           bi, alfs(bi), bts(bi)));
end
set(gca, 'FontSize', fsz)
legend(aa, 'NumColumns', 2, 'Location', 'northoutside');
xlabel('Multiplier Parameter $m$')
grid on;grid minor;grid minor;box on;
ylabel('Probability Density Function')

if savfigs
    export_fig(sprintf('./FIGS/ZC_Cdists.pdf'), '-dpdf')
end

%% * Plot
figure(1)
set(gcf, 'Color', 'white')
clf()
aa = gobjects(size(alfs));
for bi=1:length(alfs)
    aa(bi) = plot(xs, bps(:,bi), 'b-', 'LineWidth', 2); hold on
    grid on
    xlabel('Germ X')
    ylabel('PDF')
    set(gca, 'FontSize', 20, 'YTickLabel', [])
    title(sprintf('\\textbf{Beta Distribution}: $\\alpha=%.1f \\quad \\beta=%.1f$', ...
                  alfs(bi), bts(bi)));

    pause(1);
    if savfigs
        export_fig(sprintf('./FIGS/ZC_BETA_%d.eps', bi), '-depsc')
        set(aa(bi), 'Visible', 'off');
    end
    set(aa(bi), 'Visible', 'off');
end
