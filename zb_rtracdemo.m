clc
clear all
addpath('./ROUTINES/export_fig')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

savfigs = true;

%%
kt = 3;
kn = 6;
N0 = 30;
mu = 0.2;

ms = [-1 -0.25 0 0.25 1];

Nt = 32;
t = linspace(0, 2*pi, Nt+1)';
t = t(1:Nt);

uN = cos(t);
ut = sin(t);

fN = max(kn*uN+N0, 0);

nits = 3;

% Fully Stuck
fTst = mu*fN.*ms;
for mi=1:length(ms)
    for itn=1:nits
        for ti=1:Nt
            tim1 = mod(ti-1-1,Nt)+1;

            fsp = kt*(ut(ti)-ut(tim1))+fTst(tim1,mi,itn);
            if abs(fsp)<=mu*fN(ti)
                fTst(ti,mi,itn) = fsp;
            else
                fTst(ti,mi,itn) = mu*fN(ti)*sign(fsp);
            end
        end
        if itn<nits
            fTst(:,mi,itn+1) = fTst(:,mi,itn);
        end
    end
end

% Slipping
utsl = 3*sin(t);
fTsl = mu*fN.*ms;
for mi=1:length(ms)
    for itn=1:nits
        for ti=1:Nt
            tim1 = mod(ti-1-1,Nt)+1;

            fsp = kt*(utsl(ti)-utsl(tim1))+fTsl(tim1,mi,itn);
            if abs(fsp)<=mu*fN(ti)
                fTsl(ti,mi,itn) = fsp;
            else
                fTsl(ti,mi,itn) = mu*fN(ti)*sign(fsp);
            end
        end
        if itn<nits
            fTsl(:,mi,itn+1) = fTsl(:,mi,itn);
        end
    end
end


%% Plot for Slipped Case
fT = fTsl;
mpl = 1;

figure(1)
set(gcf, 'Color', 'white')
clf()
ax = axes();
plot(t, mu*fN.*[-1 1], 'k-', 'LineWidth', 2); hold on
% First Cycle
plot(t, fT(:,mpl,1), 'LineWidth', 2, 'Color', [1 1 1]*0.6)
plot(t(1:10), fT(1:10,mpl,1), 'b-', 'LineWidth', 2);
plot(t(1:2:8), fT(1:2:8,mpl,1), 'bo', 'MarkerFaceColor', 'w')
% Second Cycle
sc=plot(t, fT(:,mpl,2), 'LineWidth', 2, 'Color', 'r');

grid on
xlim([0 2*pi])
% xlim(t([1 end]))
annotation('arrow', [0.59 0.5], [0.79 0.75])
text(pi+pi/8, mu*(N0+kn), '$\mu f_N$', 'FontSize', 20)
annotation('arrow', [0.59 0.5], [0.31 0.35])
text(pi+pi/8, -mu*(N0+kn), '$-\mu f_N$', 'FontSize', 20)
set(gca, 'XTick', 0:pi/2:2*pi, ...
         'XTickLabel', {'0', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, ...
   'TickLabelInterpreter', 'latex')
xlabel('Cycle Time')
ylabel('Tangential Force (N)')
set(gca, 'FontSize', 20, 'YTickLabel', [])

pos = ax.Position;
pos(3:4) = pos(3:4) + pos(1:2);
tx = interp1(ax.XLim, pos([1 3]), t(1:10));
fy = interp1(ax.YLim, pos([2 4]), fT(1:10,mpl,1));
annotation('arrow', tx(9:10), fy(9:10), 'Color', 'b')

if savfigs
    sc.Visible = 'off';
    export_fig('./FIGS/ZB_TRDEM_0.eps', '-depsc')

    sc.Visible = 'on';
    export_fig('./FIGS/ZB_TRDEM_1.eps', '-depsc')
end


%% Plot for Sticking Case
fT = fTst;

for mpl = 1:5
    figure(1)
    set(gcf, 'Color', 'white')
    clf()
    ax = axes();
    plot(t, mu*fN.*[-1 1], 'k-', 'LineWidth', 2); hold on
    % First Cycle
    plot(t, fT(:,mpl,1), 'LineWidth', 2, 'Color', [1 1 1]*0.6)
    plot(t(1:10), fT(1:10,mpl,1), 'b-', 'LineWidth', 2);
    plot(t(1:2:8), fT(1:2:8,mpl,1), 'bo', 'MarkerFaceColor', 'w')
    % Second Cycle
    sc=plot(t, fT(:,mpl,2), 'LineWidth', 2, 'Color', 'r');

    grid on
    xlim([0 2*pi])
    % xlim(t([1 end]))
    annotation('arrow', [0.59 0.5], [0.79 0.75])
    text(pi+pi/8, mu*(N0+kn), '$\mu f_N$', 'FontSize', 20)
    annotation('arrow', [0.59 0.5], [0.31 0.35])
    text(pi+pi/8, -mu*(N0+kn), '$-\mu f_N$', 'FontSize', 20)
    set(gca, 'XTick', 0:pi/2:2*pi, ...
             'XTickLabel', {'0', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, ...
             'TickLabelInterpreter', 'latex')
    xlabel('Cycle Time')
    ylabel('Tangential Force (N)')
    set(gca, 'FontSize', 20, 'YTickLabel', [])

    pos = ax.Position;
    pos(3:4) = pos(3:4) + pos(1:2);
    tx = interp1(ax.XLim, pos([1 3]), t(1:10));
    fy = interp1(ax.YLim, pos([2 4]), fT(1:10,mpl,1));
    annotation('arrow', tx(9:10), fy(9:10), 'Color', 'b')

    if savfigs
        export_fig(sprintf('./FIGS/ZB_TRDEMst_%d.eps', mpl), '-depsc');
    end
end

%% Plot with Filled region
fT = fTst;

figure(1)
set(gcf, 'Color', 'white')
clf()
ax = axes();
plot(t, mu*fN.*[-1 1], 'k-', 'LineWidth', 2); hold on

fill(t([1:end end:-1:1]), [fT(:,1,2); fT(end:-1:1,end,2)], [1 1 1]*0.6, ...
     'FaceAlpha', 0.5);

plot(t, fT(:, [1 2 3 4 5], 2), 'r-.', 'LineWidth', 1);

grid on
xlim([0 2*pi])
% xlim(t([1 end]))
annotation('arrow', [0.59 0.5], [0.79 0.75])
text(pi+pi/8, mu*(N0+kn), '$\mu f_N$', 'FontSize', 20)
annotation('arrow', [0.59 0.5], [0.31 0.35])
text(pi+pi/8, -mu*(N0+kn), '$-\mu f_N$', 'FontSize', 20)
set(gca, 'XTick', 0:pi/2:2*pi, ...
         'XTickLabel', {'0', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, ...
         'TickLabelInterpreter', 'latex')
xlabel('Cycle Time')
ylabel('Tangential Force (N)')
set(gca, 'FontSize', 20, 'YTickLabel', [])

text(pi/2, 2.5, '$f_T\in$', 'FontSize', 20)

if savfigs
    export_fig('./FIGS/ZB_TRDEMsv.eps', '-depsc');
end
