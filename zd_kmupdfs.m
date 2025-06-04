%% * Preamble
clc
clear all
addpath('./ROUTINES/QUADRATURE/')
addpath('./ROUTINES/export_fig/');

set(0,'defaultAxesTickLabelInterpreter','default');
set(0,'defaultTextInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',13);

savfigs = true;

logsck = true;

%% * Parameters
kmu = 3e5;
ksg = kmu*.3;
if logsck
    lkmu = log(kmu^2/sqrt(kmu^2+ksg^2));
    lksg = sqrt(log(1+ksg^2/kmu^2));

    kmu = lkmu;
    ksg = lksg;
end

ktdist  = struct('mu', kmu, 'sig', ksg, 'quad', 'he_prob', ...
                 'genfun', 'ktdist.mu+ktdist.sig*xi', ...
                 'lims', [], ...
                 'pol', @(n,x) PHERM(n,x), 'pdf', 'norm');

jacfun = @(xi) exp(ktdist.mu+ktdist.sig*xi)*ktdist.sig;

Nx = 100;
xi = linspace(-4, 4, Nx);

fsz = 20;
figure(1)
poss=get(gcf, 'Position');
set(gcf, 'Color', 'white', 'Position', poss)
clf()
tiledlayout(1,1,'TileSpacing','compact', 'Padding','compact');
nexttile;

plot(exp(eval(ktdist.genfun)), pdf(ktdist.pdf,xi,0,1)./jacfun(xi), 'LineWidth', 2, ...
    'Color', 'b')
grid on;grid minor;grid minor;box on;

set(gca, 'FontSize', fsz)
xlabel('Stiffness ($\mathrm{N\,m^{-1}}$)')
ylabel('Probability Density Function')

%% * Mu
mudist = struct('ab', [0 1.25], 'mu', [], 'quad', 'legen', ...
                'genfun', 'mudist.ab(1)+(xi+1)*diff(mudist.ab)/2', ...
                'lims', [0 1.25], ...
                'pol', @(n,x) PLEGE(n,x), 'pdf', 'unif');

jacfun = @(xi) diff(mudist.ab)/2*ones(size(xi));

Nx = 100;
xi = linspace(-1, 1, Nx);

fsz = 20;
figure(2)
poss=get(gcf, 'Position');
set(gcf, 'Color', 'white', 'Position', poss)
clf()
tiledlayout(1,1,'TileSpacing','compact', 'Padding','compact');
nexttile;

plot(eval(mudist.genfun), pdf(mudist.pdf, xi, -1, 1)./jacfun(xi), 'LineWidth', 2, ...
    'Color', 'b')
set(gca, 'FontSize', fsz, 'XTick', [0:.25:1.25])

xlabel('Coefficient of Friction, $\mu$')
ylabel('Probability Density Function')
axis tight;
ylim([0 1])
grid on;grid minor;grid minor;box on;

%% * Save
if savfigs
    figure(1);
    export_fig(sprintf('./FIGS/ZD_Kpdf.pdf'), '-dpdf');

    figure(2);
    export_fig(sprintf('./FIGS/ZD_Mupdf.pdf'), '-dpdf')
end

%% * Beta distribution. Just plot
abs = [1 1; 3.5 3.5; 0.5 0.5; 5.0 2.0];

colos = DISTINGUISHABLE_COLORS(size(abs,1))*0.85;

Nx = 100;
xi = linspace(0, 1, Nx);

fsz = 20;
figure(1)
poss=get(gcf, 'Position');
set(gcf, 'Color', 'white', 'Position', poss)
clf()
tiledlayout(1,1,'TileSpacing','compact', 'Padding','compact');
nexttile;

for i=1:size(abs,1)
    plot(xi, betapdf(xi, abs(i,1), abs(i,2)), 'LineWidth', 2, 'Color', colos(i,:));
    hold on
end
grid on;grid minor;grid minor;box on;
set(gca, 'XTick', 0:0.25:1)
xlabel('$m$ Parameter')
ylabel('Probability Density Function')
set(gca, 'FontSize', fsz)
