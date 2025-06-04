%% * * Preamble
clc
clear all
addpath('./ROUTINES/')
addpath('./ROUTINES/HARMONIC/')
addpath('./ROUTINES/SOLVERS/')
addpath('./ROUTINES/QUADRATURE/')
addpath('./ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

savfigs = true;
%% * * Model Parameters (Modified, from Erhan)
pars = struct('m', 10, 'mD', 0.1, 'c', 20, ...
              'k', 1e5, 'kxy', 3e5, 'kD', 1e6);

M = diag([repmat(pars.m, 1,4) repmat(pars.mD, 1,2)]);
K = diag([repmat(pars.k, 1,4) zeros(1,2)]);

K(1:2,1:2) = K(1:2,1:2) + 3*pars.kxy*[1 -1;-1 1];
K(3:4,3:4) = K(3:4,3:4) + 1*pars.kxy*[1 -1;-1 1];
K(5:6,5:6) = K(5:6,5:6) + pars.kD*[1 -1;-1 1];
C = diag([repmat(pars.c, 1,4) zeros(1,2)]);
F0 = [0; 0; 0; 0; 0; 120]*2;
Fv = [0; -20; 5; 0; 0; 0];

Lrel1 = @(al) kron([-1 0 1], [cos(al) sin(al);-sin(al) cos(al)]);
Lrel2 = @(al) kron([0 -1 1], [cos(al) -sin(al);sin(al) cos(al)]);

ang_alpha = 40;  % Angle in degrees.
mpars = struct('M', M, 'C', C, 'K', K, ...
               'Lrels', @(al) cat(3, Lrel1(al), Lrel2(al)));

mi = 3; % Mode of interest

%% * Initialize
% 'Sim1_1111': kt,kn: Normal(3e5, 3e5*.3). mu: Unif(0, 1.25). m: beta(alf=.5, bt=.5)
% 'Sim2_1111': kt,kn: Normal(3e5, 3e5*.3). mu: Unif(0, 1.25). m: beta(alf=3.5, bt=3.5)
% 'Sim3_1111': kt,kn: Normal(3e5, 3e5*.3). mu: Unif(0, 1.25). m: beta(alf=2, bt=5)
% 'Sim4_1111': kt,kn: Normal(3e5, 3e5*.3). mu: Unif(0, 1.25). m: beta(alf=5, bt=2)
% 'Sim5_1111': kt,kn: Normal(3e5, 3e5*.3). mu: Unif(0, 1.25). m: beta(alf=.2, bt=.2)
% 'Sim6_1111': kt,kn: Normal(3e5, 3e5*.3). mu: Unif(0, 1.25). m: beta(alf=.5, bt=.2)
% 'Sim7_1111': kt,kn: Normal(3e5, 3e5*.3). mu: Unif(0, 1.25). m: beta(alf=.2, bt=.5)
% 'Sim8_1111': kt,kn: Normal(3e5, 3e5*.3). mu: Unif(0, 1.25). m: beta(alf=1, bt=1)

% Paper Cases C1 -> S8; C2 -> S2; C3 -> S1; C4 -> S4

sims = {'Sim1_1111', 'Sim2_1111', 'Sim3_1111', 'Sim4_1111', ...
        'Sim5_1111', 'Sim6_1111', 'Sim7_1111', 'Sim8_1111'};
symnam = sims{1};

% sims = {'MSim1_0011', 'MSim2_1001', 'MSim3_0101', 'MSim4_1101'};

Nqpce = 5;  % Number of quadrature points (per dimension)
for symnam=sims(4)
    % close all
    symnam = symnam{1}

    if savfigs
        figdir = ['./FIGS/' symnam '_Nq=' num2str(Nqpce)];
        mkdir(figdir)
    end    

    %% * Load Data
    uqpars = symnam(end-3:end);
    uqparsi = find(uqpars=='1');
    Nnq = sum(uqpars=='1');
    pttls = {'k_t', 'k_n', '\mu', 'm'};

    savdir = sprintf('./DATS/%s_Nq=%d', symnam, Nqpce);
    load([savdir '/sim-data.mat'], 'uqpars', 'XIs', 'WIs', ...
         'ktdist', 'kndist', 'mudist', 'mdist', 'erris', 'logsck');
    polfuns = {ktdist.pol,kndist.pol,mudist.pol,mdist.pol};
    pdffuns = {@(x) pdf(ktdist.pdf,x,0,1)./(exp(ktdist.mu+ktdist.sig*x)*ktdist.sig), ...
               @(x) pdf(kndist.pdf,x,0,1)./(exp(kndist.mu+kndist.sig*x)*kndist.sig), ...
               @(x) pdf(mudist.pdf,x,1), @(x) pdf(mdist.pdf,x,-1,1)};
    if strcmp(mudist.pdf, 'beta')
        pdffuns{3} = @(x) betapdf(x, mudist.alpha, mudist.beta);
    elseif strcmp(mudist.pdf, 'unif')
        pdffuns{3} = @(x) pdf(mudist.pdf,x,-1,1)/(diff(mudist.ab)/2);
    end
    if strcmp(mdist.pdf, 'beta')
        pdffuns{4} = @(x) betapdf(x, mdist.alpha, mdist.beta)/2;
    end

    germgens = @(xi, i) eval(ktdist.genfun).*(i==1) + ...
        eval(kndist.genfun).*(i==2) + ...
        eval(mudist.genfun).*(i==3) + ...
        eval(mdist.genfun).*(i==4);

    XIgs = zeros(Nqpce^Nnq, Nnq);
    WIgs = zeros(Nqpce^Nnq, 1);
    for isim=1:Nqpce^Nnq
        qinds = arrayfun(@(a) base2dec(a, Nqpce), dec2base(isim-1, Nqpce, Nnq))+1;

        XIgs(isim,:) = XIs(sub2ind([Nqpce Nnq], qinds, 1:Nnq));
        WIgs(isim) = prod(WIs(sub2ind([Nqpce Nnq], qinds, 1:Nnq)));
    end

    %% * Evaluate and store PCE bases
    Np = 4;  % Has to be <Nqpce since some quadrature points are actually at the roots

    IJs = cell(1,Nnq);
    [IJs{:}] = ndgrid(0:Np);
    IJs = cell2mat(cellfun(@(ij) ij(:), IJs, 'UniformOutput', false));
    IJs = IJs(sum(IJs,2)<=Np,:);
    P = size(IJs,1);
    Psi = ones(Nqpce^Nnq, P);
    Integs = ones(P, 1);  % Not really needed, just used for checking quadrature
    for n=1:P
        for i=1:Nnq
            [psi, integ] = polfuns{uqparsi(i)}(IJs(n,i), XIgs(:,i));
            Psi(:,n)  = Psi(:,n).*psi;
            Integs(n) = Integs(n)*integ;
        end
    end
    % Integs = WIgs'*(Psi.^2);
    % Integs = Integs(:);

    %% * Estimate PCE Coefficients (at each amplitude level separately)
    Ast_h = -1;
    Aen_h = 2;
    Napts = 100;  % Take 100 points equispaced in amplitude
    Qs_h = linspace(Ast_h, Aen_h, Napts)';

    LScofs = zeros(1,P);
    Wcofs = zeros(Napts,P);
    Zcofs = zeros(Napts,P);

    LSs = zeros(1, Nqpce^Nnq);
    Wss = zeros(Napts, Nqpce^Nnq);
    Zss = zeros(Napts, Nqpce^Nnq);

    lss = [];
    for isim=1:Nqpce^Nnq
        load(sprintf('%s/sim-%d.mat', savdir, isim), 'uxwL', 'xis', 'wis', 'length_scale');
        % xis should be equal to XIgs(isim,:). 
        % wis should be equal to WIgs(isim,:).

        as_h = log10((10.^uxwL(end,:))/length_scale);
        
        Ws = pchip(as_h, uxwL(end-1,:), Qs_h);
        % Zs = pchip(as_h, uxwL(end-2,:)./(2*uxwL(end-1,:)), Qs_h);
        Zs = pchip(as_h, (uxwL(end-2,:)./(2*uxwL(end-1,:))), Qs_h);

        LScofs = LScofs + WIgs(isim)*(length_scale*Psi(isim,:));
        Wcofs = Wcofs + WIgs(isim)*(Ws.*Psi(isim,:));
        Zcofs = Zcofs + WIgs(isim)*(Zs.*Psi(isim,:));

        LSs(isim) = length_scale;
        Wss(:, isim) = Ws;
        Zss(:, isim) = Zs;
    end
    LScofs = LScofs./Integs';
    Wcofs = Wcofs./Integs';
    Zcofs = Zcofs./Integs';

    Qs_h = 10.^Qs_h;

    LSvar = (LScofs(2:end).^2)*Integs(2:end);  % Length Scale variance
    Wvar = (Wcofs(:,2:end).^2)*Integs(2:end);  % Frequency variance
    Zvar = (Zcofs(:,2:end).^2)*Integs(2:end);  % Damping variance

    %% * Sobol Indices
    LS_S1 = zeros(1, Nnq);
    W_S1 = zeros(Napts, Nnq);
    Z_S1 = zeros(Napts, Nnq);
    for i=1:Nnq
        js = setdiff(1:Nnq,i);
        inds = 1+find(all(IJs(2:end,js)==0,2));  % Independent coefficients

        LS_S1(i) = (LScofs(inds).^2)*Integs(inds)./LSvar;
        W_S1(:,i) = (Wcofs(:,inds).^2)*Integs(inds)./Wvar;
        Z_S1(:,i) = (Zcofs(:,inds).^2)*Integs(inds)./Zvar;
    end

    if Nnq>=2
        INDS = find(all(IJs~=0, 2));
        inds = INDS;
        Js = nchoosek(1:Nnq,2);
        [~,jsi] = sort(Js(:,2));
        Js = Js(jsi,:);
        Ns2 = size(Js,1);
        LS_S2 = zeros(1, Ns2);
        W_S2 = zeros(Napts, Ns2);
        Z_S2 = zeros(Napts, Ns2); 
        for i=1:Ns2
            js = setdiff(1:Nnq,Js(i,:));

            inds = find(all(IJs(:,js)==0,2) & all(IJs(:,Js(i,:))~=0,2));
            LS_S2(:,i) = (LScofs(:,inds).^2)*Integs(inds)./LSvar;
            W_S2(:,i) = (Wcofs(:,inds).^2)*Integs(inds)./Wvar;
            Z_S2(:,i) = (Zcofs(:,inds).^2)*Integs(inds)./Zvar;
        end
        % W_S2 = (Wcofs(:, inds).^2)*Integs(inds)./Wvar;
        % Z_S2 = (Zcofs(:, inds).^2)*Integs(inds)./Zvar;
    end

    %% * Random Realization
    rng(1);
    Nre = 1000;
    XIre = zeros(Nre, Nnq);
    XIre = [random(ktdist.pdf, 0,1, Nre,1), random(kndist.pdf, 0,1, Nre,1), ...
            random(mudist.pdf, mudist.lims(1), mudist.lims(2), Nre,1), ...
            random(mdist.pdf, mdist.lims(1), mdist.lims(2), Nre,1)];
    if strcmp(mudist.pdf, 'beta')
        XIre(:,3) = random(mudist.pdf, mudist.alpha, mudist.beta, Nre,1);
    end
    if strcmp(mdist.pdf, 'beta')
        XIre(:,4) = random(mdist.pdf, mdist.alpha, mdist.beta, Nre,1);
    end
    XIre = XIre(:, uqparsi);
    Psire = ones(Nre, P);
    for n=1:P
        for i=1:Nnq
            [psi, integ] = polfuns{uqparsi(i)}(IJs(n,i), XIre(:,i));
            Psire(:,n)  = Psire(:,n).*psi;
        end
    end

    %% * Plot Backbone with Random Realizations + Sobol' 1
    colos = DISTINGUISHABLE_COLORS(Nnq+1);
    [~,si] = sort(vecnorm(colos'), 'ascend');
    colos = colos(si(2:end),:);
    fsz = 20;
    
    figure(1)
    poss = get(gcf, 'Position');
    set(gcf, 'Color', 'white', 'Position', [poss(1:2) 570 570])
    clf()
    tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

    recol = [[1 1 1]*0.8 0.4];
    % subplot(2,1,1)
    ax1 = nexttile;
    aa = semilogx(Qs_h, Wcofs*Psire'/2/pi, 'Color', recol); hold on
    % cc = plot(log10(Qs_h), (Wcofs(:,1)+kron([1 -1], sqrt(Wvar)))/2/pi, 'k-.', 'LineWidth', 2);
    xt = get(gca, 'XTick');
    xtl = cellfun(@(x) ['$' x '$'], get(gca, 'XTickLabels'), 'UniformOutput', false);
    bp=boxplot((Wcofs*Psire')'/2/pi, 'notch', 'on', 'positions', Qs_h, ...
               'symbol', '', 'widths', [diff(Qs_h(1:2)); diff(Qs_h)]);
    arrayfun(@(i) fill(get(bp(5,i), 'XData'), get(bp(5,i), 'YData'), ...
                       'w', 'EdgeColor', 'b'), 1:size(bp,2));
    arrayfun(@(i) plot(get(bp(6,i), 'XData'), get(bp(6,i), 'YData'), 'r-'), 1:size(bp,2));
    % set(bp(3:end,:), 'LineWidth', 1.5)
    bb = plot(Qs_h, Wcofs(:,1)/2/pi, 'k-', 'LineWidth', 2);
    grid on; grid minor; grid minor
    set(gca, 'XTick', xt, 'XTickLabels', xtl, 'TickLabelinterpreter', 'latex', ...
             'XScale', 'log');
    xlim(10.^[Ast_h Aen_h])
    ylabel('Frequency [Hz]')
    set(gca, 'FontSize', fsz);

    zsc = 100;
    % zsc = 1.0;
    % subplot(2,1,2)
    ax2 = nexttile;
    semilogx(Qs_h, Zcofs*Psire'*zsc, 'Color', recol); hold on
    % plot(log10(Qs_h), (Zcofs(:,1)+kron([1 -1], sqrt(Zvar)))*100, 'k-.', 'LineWidth', 2)
    xt = get(gca, 'XTick');
    xtl = cellfun(@(x) ['$' x '$'], get(gca, 'XTickLabels'), 'UniformOutput', false);
    bp=boxplot((Zcofs*Psire')'*zsc, 'notch', 'on', 'positions', Qs_h, ...
               'symbol', '', 'widths', [diff(Qs_h(1:2)); diff(Qs_h)]);
    arrayfun(@(i) fill(get(bp(5,i), 'XData'), get(bp(5,i), 'YData'), ...
                       'w', 'EdgeColor', 'b'), 1:size(bp,2));
    arrayfun(@(i) plot(get(bp(6,i), 'XData'), get(bp(6,i), 'YData'), 'r-'), 1:size(bp,2));
    % set(bp(3:end,:), 'LineWidth', 1.5);
    semilogx(Qs_h, Zcofs(:,1)*zsc, 'k-', 'LineWidth', 2)
    set(gca, 'XTick', xt, 'XTickLabels', xtl, 'TickLabelinterpreter', 'latex', ...
             'XScale', 'log');
    grid on; grid minor; grid minor
    ylabel('Damping Factor [\%]')
    xlabel('Scaled Modal Amplitude $q/\ell$')
    xlim(10.^[Ast_h Aen_h])
    yl = ylim;
    ylim([-1.5 max(yl)])
    set(gca, 'FontSize', fsz);

    if savfigs
    end

    axes(ax1)
    yyaxis right;
    plot(Qs_h, W_S1(:,1), '-', 'LineWidth', 2, 'Color', colos(1,:))
    if Nnq>=2
        plot(Qs_h, W_S1(:,2), '-', 'LineWidth', 2, 'Color', colos(2,:))
    end
    switch Nnq
      case 1
      case 2
        plot(Qs_h, W_S2, '-', 'LineWidth', 2, 'Color', colos(3,:))
      case 3
        plot(Qs_h, W_S1(:,3), '-', 'LineWidth', 2, 'Color', colos(3,:))
      case 4
        plot(Qs_h, W_S1(:,3), '-', 'LineWidth', 2, 'Color', colos(3,:))
        plot(Qs_h, W_S1(:,4), '-', 'LineWidth', 2, 'Color', colos(4,:))
      otherwise
        warning('Nnq>3')
    end
    ylabel('Sobol'' Index')
    set(gca, 'YTick', 0:0.2:1)
    ylim([-0.01 1.01])

    axes(ax2)
    yyaxis right;
    aa = plot(Qs_h, Z_S1(:,1), '-', 'LineWidth', 2, 'Color', colos(1,:));
    if Nnq>=2
        bb = plot(Qs_h, Z_S1(:,2), '-', 'LineWidth', 2, 'Color', colos(2,:));
    end
    switch Nnq
      case 1
      case 2
        cc = plot(Qs_h, Z_S2, '-', 'LineWidth', 2, 'Color', colos(3,:));
        legend([aa bb cc], sprintf('$SU_{%s}$',pttls{uqparsi(1)}), ...
               sprintf('$SU_{%s}$',pttls{uqparsi(2)}), ...
               sprintf('$SU_{%s,%s}$',pttls{uqparsi(1)},pttls{uqparsi(2)}), ...
               'Location', 'northeast')
      case 3
        cc = plot(Qs_h, Z_S1(:,3), '-', 'LineWidth', 2);
        legend([aa bb cc], sprintf('$SU_{%s}$',pttls{uqparsi(1)}), ...
               sprintf('$SU_{%s}$',pttls{uqparsi(2)}), ...
               sprintf('$SU_{%s}$',pttls{uqparsi(3)}), ...
               'Location', 'northeast')
      case 4
        cc = plot(Qs_h, Z_S1(:,3), '-', 'LineWidth', 2, 'Color', colos(3,:));
        dd = plot(Qs_h, Z_S1(:,4), '-', 'LineWidth', 2, 'Color', colos(4,:));
        ll = legend([aa bb cc dd], sprintf('$SU_{%s}$',pttls{uqparsi(1)}), ...
                    sprintf('$SU_{%s}$',pttls{uqparsi(2)}), ...
                    sprintf('$SU_{%s}$',pttls{uqparsi(3)}), ...
                    sprintf('$SU_{%s}$',pttls{uqparsi(4)}), ...
                    'Location', 'best');
      otherwise
        warning('Nnq>3');
    end
    ylabel('Sobol'' Index')
    set(gca, 'YTick', 0:0.2:1)
    ylim([-0.01 1.01])

    if savfigs
        % export_fig([figdir '/D2_BBwS1.png'], '-dpng', '-r300');
        export_fig([figdir '/D2_BBwS1.png'], '-dpng');
    end

    %% * Plot Backbone with Random Realizations + Sobol' 2
    if Nnq==4
        colos = DISTINGUISHABLE_COLORS(3);

        figure(10)
        poss = get(gcf, 'Position');
        set(gcf, 'Color', 'white', 'Position', [poss(1:2) 570 570])
        clf()
        tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

        recol = [[1 1 1]*0.8 0.4];
        ax1 = nexttile;
        aa = semilogx(Qs_h, Wcofs*Psire'/2/pi, 'Color', recol); hold on
        xt = get(gca, 'XTick');
        xtl = cellfun(@(x) ['$' x '$'], get(gca, 'XTickLabels'), 'UniformOutput', false);
        bp=boxplot((Wcofs*Psire')'/2/pi, 'notch', 'on', 'positions', Qs_h, ...
                   'symbol', '', 'widths', [diff(Qs_h(1:2)); diff(Qs_h)]);
        arrayfun(@(i) fill(get(bp(5,i), 'XData'), get(bp(5,i), 'YData'), ...
                           'w', 'EdgeColor', 'b'), 1:size(bp,2));
        arrayfun(@(i) plot(get(bp(6,i), 'XData'), get(bp(6,i), 'YData'), 'r-'), 1:size(bp,2));
        bb = plot(Qs_h, Wcofs(:,1)/2/pi, 'k-', 'LineWidth', 2);
        grid on; grid minor; grid minor
        set(gca, 'XTick', xt, 'XTickLabels', xtl, 'TickLabelinterpreter', 'latex', ...
                 'XScale', 'log');
        xlim(10.^[Ast_h Aen_h])
        ylabel('Frequency [Hz]')
        set(gca, 'FontSize', fsz);

        ax2 = nexttile;
        semilogx(Qs_h, Zcofs*Psire'*zsc, 'Color', recol); hold on
        xt = get(gca, 'XTick');
        xtl = cellfun(@(x) ['$' x '$'], get(gca, 'XTickLabels'), 'UniformOutput', false);
        bp=boxplot((Zcofs*Psire')'*zsc, 'notch', 'on', 'positions', Qs_h, ...
                   'symbol', '', 'widths', [diff(Qs_h(1:2)); diff(Qs_h)]);
        arrayfun(@(i) fill(get(bp(5,i), 'XData'), get(bp(5,i), 'YData'), ...
                           'w', 'EdgeColor', 'b'), 1:size(bp,2));
        arrayfun(@(i) plot(get(bp(6,i), 'XData'), get(bp(6,i), 'YData'), 'r-'), 1:size(bp,2));
        semilogx(Qs_h, Zcofs(:,1)*zsc, 'k-', 'LineWidth', 2)
        set(gca, 'XTick', xt, 'XTickLabels', xtl, 'TickLabelinterpreter', 'latex', ...
                 'XScale', 'log');
        grid on; grid minor; grid minor
        ylabel('Damping Factor [\%]')
        xlabel('Scaled Modal Amplitude $q/\ell$')
        xlim(10.^[Ast_h Aen_h])
        yl = ylim;
        ylim([-1.5 max(yl)])

        set(gca, 'FontSize', fsz);

        for ipi=1:2
            if (ipi>1)
                delete([qa qb qc aa bb cc ll]);
            end
            
            axes(ax1)
            yyaxis right;
            qa = plot(Qs_h, W_S2(:,(ipi-1)*3+1), '-', 'LineWidth', 2, ...
                      'Color', colos(1,:));
            qb = plot(Qs_h, W_S2(:,(ipi-1)*3+2), '-', 'LineWidth', 2, ...
                      'Color', colos(2,:));
            qc = plot(Qs_h, W_S2(:,(ipi-1)*3+3), '-', 'LineWidth', 2, ...
                      'Color', colos(3,:));
            ylabel('Sobol'' Index')
            % set(gca, 'YTick', 0:0.2:1)
            % ylim([-0.01 1.01])

            axes(ax2)
            yyaxis right;
            aa = plot(Qs_h, Z_S2(:,(ipi-1)*3+1), '-', 'LineWidth', 2, ...
                      'Color', colos(1,:));
            bb = plot(Qs_h, Z_S2(:,(ipi-1)*3+2), '-', 'LineWidth', 2, ...
                      'Color', colos(2,:));
            cc = plot(Qs_h, Z_S2(:,(ipi-1)*3+3), '-', 'LineWidth', 2, ...
                      'Color', colos(3,:));
            ylabel('Sobol'' Index')
            ll = legend([aa bb cc], ...
                        sprintf('$SU_{%s,%s}$', ...
                                pttls{uqparsi(Js((ipi-1)*3+1,1))}, ...
                                pttls{uqparsi(Js((ipi-1)*3+1,2))}), ...
                        sprintf('$SU_{%s,%s}$', ...
                                pttls{uqparsi(Js((ipi-1)*3+2,1))}, ...
                                pttls{uqparsi(Js((ipi-1)*3+2,2))}), ...
                        sprintf('$SU_{%s,%s}$', ...
                                pttls{uqparsi(Js((ipi-1)*3+3,1))}, ...
                                pttls{uqparsi(Js((ipi-1)*3+3,2))}), ...
                        'Location', 'northwest');

            if savfigs
                figdir = ['./FIGS/' symnam '_Nq=' num2str(Nqpce)];
                mkdir(figdir)

                export_fig([figdir '/D2_BBwrS2_' num2str(ipi) '.png'], '-dpng', '-r300');
            end
        end        
    end

    %% * Plot Length_Scale
    figure(2)
    poss=get(gcf, 'Position');
    set(gcf, 'Color', 'white', 'Position', poss)
    clf()
    tiledlayout(1,1,'TileSpacing','compact', 'Padding','compact');
    nexttile;
    wd = 0.075;

    aa = plot(1.05*ones(size(LSs)), LSs, '.', 'Color', [1 1 1]*0.6, 'MarkerSize', 15); hold on
    % bb = plot(1.05+[-1 1]*wd, LScofs(1)*[1 1], 'r-', 'LineWidth', 2);
    % cc = plot(1.05+[-1 1 1 -1 -1]*wd, LScofs(1)+[-1 -1 1 1 -1]*sqrt(LSvar), 'b-', 'LineWidth', 2);
    bb=boxplot(LSs, 'positions', 1.05);
    set(bb, 'LineWidth', 2);
    xlim([0.8 1.8])
    ylabel('Length Scale $\ell$ [$\mathrm{kg^{1/2}\,m}$]')
    % set(gca, 'YScale', 'log')
    grid on; grid minor; grid minor

    yyaxis right
    xxss = linspace(mean(xlim), max(xlim), Nnq+2);
    stem(xxss(2:end-1), LS_S1, 'LineWidth', 2, 'MarkerFaceColor', 'w'); hold on
    set(gca, 'XTick', xxss(2:end-1), 'XTickLabels', ...
             cellfun(@(pt) sprintf('$%s$', pt), pttls(uqparsi), ...
                     'UniformOutput', false), 'TickLabelInterpreter', 'latex')
    ylabel('Sobol'' Index')
    ylim([0 1]);
    plot(mean(xlim)*[1 1], ylim, 'k-')
    grid on; grid minor; grid minor
    set(gca, 'FontSize', fsz);

    if savfigs
        % export_fig([figdir '/D2_LS.eps'], '-depsc');

        export_fig([figdir '/D2_LS.png'], '-dpng', '-r300');
    end

    %% * Plotting LS with S2
    figure(20)
    poss=get(gcf, 'Position');
    set(gcf, 'Color', 'white', 'Position', [poss(1:2) 700 400])
    clf()
    tiledlayout(1,1,'TileSpacing','compact', 'Padding','compact');
    nexttile;

    wd = 0.025;

    aa = plot(1.05*ones(size(LSs)), LSs, '.', 'Color', [1 1 1]*0.6, 'MarkerSize', 15); hold on
    % bb = plot(1.05+[-1 1]*wd, LScofs(1)*[1 1], 'r-', 'LineWidth', 2);
    % cc = plot(1.05+[-1 1 1 -1 -1]*wd, LScofs(1)+[-1 -1 1 1 -1]*sqrt(LSvar), 'b-', 'LineWidth', 2);

    bb=boxplot(LSs, 'positions', 1.05);
    set(bb, 'LineWidth', 2);
    xlim([0.95 1.8])
    ylabel('Length Scale (m)')
    % set(gca, 'YScale', 'log')
    grid on; grid minor; grid minor;

    yyaxis right
    if Nnq>2
        xxss = linspace(1.15, max(xlim), nchoosek(Nnq,2)+2);
        stem(xxss(2:end-1), LS_S2, 'LineWidth', 2, 'MarkerFaceColor', 'w'); hold on
        set(gca, 'XTick', xxss(2:end-1), 'XTickLabels', ...
                 cellfun(@(pt1,pt2) sprintf('$(%s,%s)$', pt1,pt2), ...
                         pttls(uqparsi(Js(:,1))), pttls(uqparsi(Js(:,2))), ...
                         'UniformOutput', false), 'TickLabelInterpreter', 'latex')
        ylabel('Sobol'' Index')
        plot(xxss(1)*[1 1], ylim, 'k-')
        grid on; grid minor; grid minor;
        set(gca, 'FontSize', 16)
    end

    if savfigs
        export_fig([figdir '/D2_LS2.png'], '-dpng', '-r300');
    end


    %% * Plotting PDFs
    sc = 1.5;
    uns = {'[$\mathrm{kN\,m^{-1}}$]', '[$\mathrm{kN\,m^{-1}}$]', '', ''};

    fsz = 20;
    
    i = 2;
    for i=1:Nnq
        figure(i*10)
        poss=get(gcf, 'Position');
        set(gcf, 'Color', 'white', 'Position', [poss(1:2) 560 470])
        clf()
        tiledlayout(1,1,'TileSpacing','compact', 'Padding','compact');
        nexttile;

        mnmx = [-inf inf];
        if uqparsi(i)==4
            mnmx = [-1 1];
        end
        if uqparsi(i)==3 && strcmp(mudist.pdf,'unif')
            mnmx = [-1 1];
        end
        if min(XIs(:,i))>0
            Xs = linspace(max(min(XIs(:,i))/sc,mnmx(1)), ...
                          min(max(XIs(:,i))*sc,mnmx(2)), 100);
        else
            Xs = linspace(max(min(XIs(:,i))*sc,mnmx(1)), ...
                          min(max(XIs(:,i))*sc,mnmx(2)), 100);
        end
        if uqparsi(i)==4
            Xs = linspace(0, 1, 100);
        end        
        if (uqparsi(i)==1 || uqparsi(i)==2) && logsck
            plot(exp(germgens(Xs,uqparsi(i))), pdffuns{uqparsi(i)}(Xs), ...
                 'LineWidth', 2, 'Color', 'b'); hold on
            stem(exp(germgens(XIs(:,i),uqparsi(i))), pdffuns{uqparsi(i)}(XIs(:,i)), ...
                 'MarkerFaceColor', 'w', 'LineWidth', 1.5, 'Color', 'r')
            % set(gca, 'XScale', 'log')
        else
            plot(germgens(Xs,uqparsi(i)), pdffuns{uqparsi(i)}(Xs), ...
                 'LineWidth', 2, 'Color', 'b'); hold on
            stem(germgens(XIs(:,i),uqparsi(i)), pdffuns{uqparsi(i)}(XIs(:,i)), ...
                 'MarkerFaceColor', 'w', 'LineWidth', 1.5, 'Color', 'r')
        end
        % xlabel(sprintf('Germ $x_{%s}$', pttls{uqparsi(i)}))
        xlabel(sprintf('Factor Value $%s$ %s', pttls{uqparsi(i)}, uns{uqparsi(i)}))
        ylabel('Probability Density Function')

        switch uqparsi(i)
          case 3
            xlim([0 1.25]);
            set(gca, 'XTick', 0:0.25:1.25)
          case 4
            xlim([-1 1]);
        end
        

        grid on; grid minor; grid minor;
        set(gca, 'FontSize', fsz)

        if savfigs
            export_fig([figdir '/D2_PDF_' num2str(i) '.pdf'], '-dpdf');
        end
    end

    %% * If Nnq is 2, plot the two germs
    fac_gIs = cell2mat(arrayfun(@(i) germgens(XIgs(:,i),uqparsi(i)), 1:Nnq, 'UniformOutput', false));
    if logsck
        for i=1:Nnq
            if (uqparsi(i)==1 || uqparsi(i)==2)
                fac_gIs(:, i) = exp(fac_gIs(:, i));
            end
        end
    end
    if Nnq==2
        uns = {'kN/m', 'kN/m', '', ''};
        
        figure(100)
        set(gcf, 'Color', 'white')
        clf()
        plot(fac_gIs(:,1), fac_gIs(:,2), 'ko', 'MarkerFaceColor', 'k');
        xlabel(sprintf('Factor Value $%s$ %s', pttls{uqparsi(1)}, uns{uqparsi(1)}))
        ylabel(sprintf('Factor Value $%s$ %s', pttls{uqparsi(2)}, uns{uqparsi(2)}))
        grid on

        if savfigs
            export_fig([figdir '/D2_PGRID.pdf'], '-dpdf');
        end
    end

    %% * Plot all the realizations used for the quadrature
    % load([savdir '/sim-data.mat'], 'erris');

    % colos = DISTINGUISHABLE_COLORS(Nqpce^(Nnq-1));
    % scs = [1e-3 1e-3 1.0 1.0];
    % uns = {'kN/m', 'kN/m', '', ''};

    % figure(101)
    % clf()
    % set(gcf, 'Color', 'white')
    % aa = gobjects(Nqpce^(Nnq-1),1);
    % for isim=find(~erris)'
    %     load(sprintf('%s/sim-%d.mat', savdir, isim), 'uxwL', 'xis', 'wis', 'length_scale');

    %     subplot(2,1,1)
    %     semilogx(10.^uxwL(end,:)/length_scale, ...
    %              uxwL(end-1,:)/2/pi, '-', 'Color', ...
    %              colos(floor((isim-1)/Nqpce)+1,:)); hold on

    %     subplot(2,1,2)
    %     aa(floor((isim-1)/Nqpce)+1) = semilogx(10.^uxwL(end,:)/length_scale, ...
    %                                            uxwL(end-2,:)./(2*uxwL(end-1,:))*100, ...
    %                                            '-', 'Color', ...
    %                                            colos(floor((isim-1)/Nqpce)+1,:)); hold on
    %     switch Nnq
    %       case 1
    %       case 2
    %         ll = legend(aa(floor((isim-1)/Nqpce)+1), ...
    %                     sprintf('$%s$ = %.2f %s', pttls{uqparsi(1)}, ...
    %                             fac_gIs(isim,1)*scs(uqparsi(1)), uns{uqparsi(1)}));
    %       case 3
    %         ll = legend(aa(floor((isim-1)/Nqpce)+1), ...
    %                     sprintf('$(%s, %s)$ = (%.2f %s, %.2f %s) ', ...
    %                             pttls{uqparsi(1)}, pttls{uqparsi(2)}, ...
    %                             fac_gIs(isim,1)*scs(uqparsi(1)), uns{uqparsi(1)}, ...
    %                             fac_gIs(isim,2)*scs(uqparsi(2)), uns{uqparsi(2)}));
    %       otherwise
    %         warning('Nnq>3')
    %     end    
    %     % drawnow
    % end
    % subplot(2,1,1)
    % xlim(10.^[Ast_h Aen_h])
    % if Nnq~=4
    %     ll=legend(aa, 'Location', 'northwest');
    % end
    % ylabel('Natural Frequency (Hz)')
    % subplot(2,1,2)
    % xlim(10.^[Ast_h Aen_h])
    % ylabel('Damping Factor (\%)')
    % xlabel('Modal Amplitude')

    % % if savfigs
    % %     if Nqpce^(Nnq-1)>6
    % %         set(ll, 'Visible', 'off');
    % %     end
    % %     export_fig([figdir '/D2_BB_data.eps'], '-depsc');
    % % end

    %% * Plotting PCE Results
    % figure(102)
    % poss = get(gcf, 'Position');
    % set(gcf, 'Color', 'white', 'Position', [poss(1:2) 570 500])
    % clf()
    % tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');
    % ax1 = nexttile;
    % ax2 = nexttile;
    % for isim=1:Nqpce^Nnq
    %     load(sprintf('%s/sim-%d.mat', savdir, isim), ...
    %          'uxwL', 'xis', 'wis', 'length_scale');

    %     % subplot(2,1,1)
    %     axes(ax1)
    %     aa = semilogx(10.^uxwL(end,:)/length_scale, uxwL(end-1,:)/2/pi, '-', ...
    %                   'Color', [1 1 1]*0.6, 'LineWidth', 2); hold on

    %     % subplot(2,1,2)
    %     axes(ax2)
    %     semilogx(10.^uxwL(end,:)/length_scale, uxwL(end-2,:)./(2*uxwL(end-1,:))*100, '-', ...
    %              'Color', [1 1 1]*0.6, 'LineWidth', 2); hold on
    %     % drawnow
    % end
    % % subplot(2,1,1)
    % axes(ax1);
    % bb = plot(Qs_h, Wcofs(:,1)/2/pi, 'k-', 'LineWidth', 2);
    % cc = plot(Qs_h, (Wcofs(:,1)+kron([1 -1], sqrt(Wvar)))/2/pi, 'k-.', 'LineWidth', 2);

    % ylabel('Natural Frequency (Hz)')

    % yyaxis right;
    % plot(Qs_h, W_S1(:,1), 'b-', 'LineWidth', 2)
    % if Nnq>=2
    %     plot(Qs_h, W_S1(:,2), 'r-', 'LineWidth', 2)
    % end
    % switch Nnq
    %   case 1
    %   case 2
    %     plot(Qs_h, W_S2, 'm-', 'LineWidth', 2)
    %   case 3
    %     plot(Qs_h, W_S1(:,3), 'm-', 'LineWidth', 2)
    %   case 4
    %     plot(Qs_h, W_S1(:,3), 'm-', 'LineWidth', 2)
    %     plot(Qs_h, W_S1(:,4), 'c-', 'LineWidth', 2)
    %   otherwise
    %     warning('Nnq>3')
    % end

    % ylabel('Sobol'' Index')
    % ylim([-0.01 1.01])

    % yyaxis left
    % legend([aa bb cc(1)], 'Simulations', 'PCE Mean', 'PCE $\mu\pm\sigma$', 'Location', 'northeast')

    % xlim(10.^[Ast_h Aen_h])
    % % subplot(2,1,2)
    % axes(ax2);
    % plot(Qs_h, Zcofs(:,1)*100, 'k-', 'LineWidth', 2)
    % plot(Qs_h, (Zcofs(:,1)+kron([1 -1], sqrt(Zvar)))*100, 'k-.', 'LineWidth', 2)
    % ylabel('Damping Factor (\%)')
    % xlabel('Modal Amplitude')
    % xlim(10.^[Ast_h Aen_h])

    % yyaxis right;
    % aa = plot(Qs_h, Z_S1(:,1), 'b-', 'LineWidth', 2);
    % if Nnq>=2
    %     bb = plot(Qs_h, Z_S1(:,2), 'r-', 'LineWidth', 2);
    % end
    % switch Nnq
    %   case 1
    %   case 2
    %     cc = plot(Qs_h, Z_S2, 'm-', 'LineWidth', 2);
    %     legend([aa bb cc], sprintf('$SU_{%s}$',pttls{uqparsi(1)}), ...
    %            sprintf('$SU_{%s}$',pttls{uqparsi(2)}), ...
    %            sprintf('$SU_{%s,%s}$',pttls{uqparsi(1)},pttls{uqparsi(2)}), ...
    %            'Location', 'northeast')
    %   case 3
    %     cc = plot(Qs_h, Z_S1(:,3), 'm-', 'LineWidth', 2);
    %     legend([aa bb cc], sprintf('$SU_{%s}$',pttls{uqparsi(1)}), ...
    %            sprintf('$SU_{%s}$',pttls{uqparsi(2)}), ...
    %            sprintf('$SU_{%s}$',pttls{uqparsi(3)}), ...
    %            'Location', 'northeast')
    %   case 4
    %     cc = plot(Qs_h, Z_S1(:,3), 'm-', 'LineWidth', 2);
    %     dd = plot(Qs_h, Z_S1(:,4), 'c-', 'LineWidth', 2);
    %     legend([aa bb cc dd], sprintf('$SU_{%s}$',pttls{uqparsi(1)}), ...
    %            sprintf('$SU_{%s}$',pttls{uqparsi(2)}), ...
    %            sprintf('$SU_{%s}$',pttls{uqparsi(3)}), ...
    %            sprintf('$SU_{%s}$',pttls{uqparsi(4)}), ...
    %            'Location', 'northoutside', 'numcolumns', 4)
    %   otherwise
    %     warning('Nnq>3');
    % end
    % ylabel('Sobol'' Index')
    % ylim([-0.01 1.01])

    % if savfigs
    %     figdir = ['./FIGS/' symnam '_Nq=' num2str(Nqpce)];
    %     mkdir(figdir)
    %     % export_fig([figdir '/D2_BB.eps'], '-depsc');

    %     export_fig([figdir '/D2_BB.png'], '-dpng', '-r300');
    % end
    %% * 
end
