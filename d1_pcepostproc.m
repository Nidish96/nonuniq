clc
% close all
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

savfig = false;
% * Model Parameters (Modified, from Erhan)
pars = struct('m', 10, 'mD', 0.1, 'c', 20, ...
    'k', 1e5, 'kxy', 3e5, 'kD', 1e6);

M = diag([repmat(pars.m, 1,4) repmat(pars.mD, 1,2)]);
K = diag([repmat(pars.k, 1,4) zeros(1,2)]);

K(1:2,1:2) = K(1:2,1:2) + 3*pars.kxy*[1 -1;-1 1];
K(3:4,3:4) = K(3:4,3:4) + 1*pars.kxy*[1 -1;-1 1];
K(5:6,5:6) = K(5:6,5:6) + pars.kD*[1 -1;-1 1];
C = diag([repmat(pars.c, 1,4) zeros(1,2)]);
F0 = [0; 0; 0; 0; 0; 120];
Fv = [0; -20; 5; 0; 0; 0];

Lrel1 = @(al) kron([-1 0 1], [cos(al) sin(al);-sin(al) cos(al)]);
Lrel2 = @(al) kron([0 -1 1], [cos(al) -sin(al);sin(al) cos(al)]);

ang_alpha = 40;  % Angle in degrees.
mpars = struct('M', M, 'C', C, 'K', K, ...
    'Lrels', @(al) cat(3, Lrel1(al), Lrel2(al)));

mi = 3; % Mode of interest

% * Initialize 
sims = {'MSim1_0011', 'MSim2_1001', 'MSim3_0101', 'MSim4_1101', 'Bsim4_0011'};
symnam = sims{5};
Nqpce = 5;  % Number of quadrature points (per dimension)

uqpars = symnam(end-3:end);
uqparsi = find(uqpars=='1');
Nnq = sum(uqpars=='1');
pttls = {'k_t', 'k_n', '\mu', 'm'};

savdir = sprintf('./DATS/%s_Nq=%d', symnam, Nqpce);
load([savdir '/sim-data.mat'], 'uqpars', 'XIs', 'WIs', 'ktdist', 'kndist', 'mudist', 'mdist', 'erris');
polfuns = {ktdist.pol,kndist.pol,mudist.pol,mdist.pol};
pdffuns = {@(x) pdf(ktdist.pdf,x,0,1), @(x) pdf(kndist.pdf,x,0,1), @(x) pdf(mudist.pdf,x,1), @(x) pdf(mdist.pdf,x,-1,1)};

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

% * Evaluate and store PCE bases
Np = 4;  % Has to be lesser than Nqpce since some quadrature points are actually at the roots

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

% * Estimate PCE Coefficients (at each amplitude level separately)
Ast = -5;
Aen = -1.5;
Napts = 100;  % Take 100 points equispaced in amplitude
Qs = linspace(Ast, Aen, Napts)';

Wcofs = zeros(Napts,P);
Zcofs = zeros(Napts,P);

Wss = zeros(Napts, Nqpce^Nnq);
Zss = zeros(Napts, Nqpce^Nnq);
for isim=1:Nqpce^Nnq
    load(sprintf('%s/sim-%d.mat', savdir, isim), 'uxwL', 'xis', 'wis');
    % xis should be equal to XIgs(isim,:). 
    % wis should be equal to WIgs(isim,:).
    
    Ws = pchip(uxwL(end,:),uxwL(end-1,:), Qs);
    Zs = pchip(uxwL(end,:),uxwL(end-2,:)./(2*uxwL(end-1,:)), Qs);

    Wcofs = Wcofs + WIgs(isim)*(Ws.*Psi(isim,:));
    Zcofs = Zcofs + WIgs(isim)*(Zs.*Psi(isim,:));

    Wss(:, isim) = Ws;
    Zss(:, isim) = Zs;
end
Wcofs = Wcofs./Integs';
Zcofs = Zcofs./Integs';

Qs = 10.^Qs;

Wvar = (Wcofs(:,2:end).^2)*Integs(2:end);  % Frequency variance
Zvar = (Zcofs(:,2:end).^2)*Integs(2:end);  % Damping variance

% * Sobol Indices
W_S1 = zeros(Napts, Nnq);
Z_S1 = zeros(Napts, Nnq);
for i=1:Nnq
    js = setdiff(1:Nnq,i);
    inds = 1+find(all(IJs(2:end,js)==0,2));  % Independent coefficients

    W_S1(:,i) = (Wcofs(:,inds).^2)*Integs(inds)./Wvar;
    Z_S1(:,i) = (Zcofs(:,inds).^2)*Integs(inds)./Zvar;
end

INDS = find(all(IJs~=0, 2));
inds = INDS;
Js = nchoosek(1:Nnq,2);
Ns2 = size(Js,1);
W_S2 = zeros(Napts, Ns2);
Z_S2 = zeros(Napts, Ns2); 
for i=1:Ns2
    js = setdiff(1:Nnq,Js(i,:));

    inds = find(all(IJs(:,js)==0,2) & all(IJs(:,Js(i,:))~=0,2));
    W_S2(:,i) = (Wcofs(:,inds).^2)*Integs(inds)./Wvar;
    Z_S2(:,i) = (Zcofs(:,inds).^2)*Integs(inds)./Zvar;
end
% W_S2 = (Wcofs(:, inds).^2)*Integs(inds)./Wvar;
% Z_S2 = (Zcofs(:, inds).^2)*Integs(inds)./Zvar;

% * Plotting PCE Results
figure(1)
clf()
set(gcf, 'Color', 'white')
for isim=1:Nqpce^Nnq
    load(sprintf('%s/sim-%d.mat', savdir, isim), 'uxwL', 'xis', 'wis');

    subplot(2,1,1)
    aa = semilogx(10.^uxwL(end,:), uxwL(end-1,:)/2/pi, '-', ...
                  'Color', [1 1 1]*0.6, 'LineWidth', 2); hold on
    
    subplot(2,1,2)
    semilogx(10.^uxwL(end,:), uxwL(end-2,:)./(2*uxwL(end-1,:))*100, '-', ...
             'Color', [1 1 1]*0.6, 'LineWidth', 2); hold on
    drawnow
end
subplot(2,1,1)
bb = plot(Qs, Wcofs(:,1)/2/pi, 'k-', 'LineWidth', 2);
cc = plot(Qs, (Wcofs(:,1)+kron([1 -1], sqrt(Wvar)))/2/pi, 'k-.', 'LineWidth', 2);
ylabel('Natural Frequency (Hz)')

yyaxis right;
plot(Qs, W_S1(:,1), 'b-', 'LineWidth', 2)
plot(Qs, W_S1(:,2), 'r-', 'LineWidth', 2)
switch Nnq
    case 2
        plot(Qs, W_S2, 'm-', 'LineWidth', 2)
    case 3
        plot(Qs, W_S1(:,3), 'm-', 'LineWidth', 2)
    otherwise
        warning('Nnq>3')
end

ylabel('Sobol'' Index')
ylim([-0.01 1.01])

yyaxis left
legend([aa bb cc(1)], 'Simulations', 'PCE Mean', 'PCE $\mu\pm\sigma$', 'Location', 'northeast')

xlim(10.^[Ast Aen])
subplot(2,1,2)
plot(Qs, Zcofs(:,1)*100, 'k-', 'LineWidth', 2)
plot(Qs, (Zcofs(:,1)+kron([1 -1], sqrt(Zvar)))*100, 'k-.', 'LineWidth', 2)
ylabel('Damping Factor (\%)')
xlabel('Modal Amplitude')
xlim(10.^[Ast Aen])

yyaxis right;
aa = plot(Qs, Z_S1(:,1), 'b-', 'LineWidth', 2);
bb = plot(Qs, Z_S1(:,2), 'r-', 'LineWidth', 2);
switch Nnq
    case 2
        cc = plot(Qs, Z_S2, 'm-', 'LineWidth', 2);
        legend([aa bb cc], sprintf('$SU_{%s}$',pttls{uqparsi(1)}), ...
            sprintf('$SU_{%s}$',pttls{uqparsi(2)}), ...
            sprintf('$SU_{%s,%s}$',pttls{uqparsi(1)},pttls{uqparsi(2)}), ...
            'Location', 'northeast')
    case 3
        cc = plot(Qs, Z_S1(:,3), 'm-', 'LineWidth', 2);
        legend([aa bb cc], sprintf('$SU_{%s}$',pttls{uqparsi(1)}), ...
            sprintf('$SU_{%s}$',pttls{uqparsi(2)}), ...
            sprintf('$SU_{%s}$',pttls{uqparsi(3)}), ...
            'Location', 'northeast')
    otherwise
        warning('Nnq>3');
end
ylabel('Sobol'' Index')
ylim([-0.01 1.01])

if savfig
    figdir = ['./FIGS/' symnam '_Nq=' num2str(Nqpce)];
    mkdir(figdir)
    export_fig([figdir '/BB.eps'], '-depsc');

    export_fig([figdir '/BB.png'], '-dpng', '-r300');
end

% *
if Nnq==3
    figure(2)
    clf()
    set(gcf, 'Color', 'white')
    for isim=1:Nqpce^Nnq
        load(sprintf('%s/sim-%d.mat', savdir, isim), 'uxwL', 'xis', 'wis');

        subplot(2,1,1)
        aa = semilogx(10.^uxwL(end,:), uxwL(end-1,:)/2/pi, '-', ...
                      'Color', [1 1 1]*0.6, 'LineWidth', 2); hold on
        
        subplot(2,1,2)
        semilogx(10.^uxwL(end,:), uxwL(end-2,:)./(2*uxwL(end-1,:))*100, '-', ...
                 'Color', [1 1 1]*0.6, 'LineWidth', 2); hold on
        drawnow
    end
    subplot(2,1,1)
    bb = plot(Qs, Wcofs(:,1)/2/pi, 'k-', 'LineWidth', 2);
    cc = plot(Qs, (Wcofs(:,1)+kron([1 -1], sqrt(Wvar)))/2/pi, 'k-.', 'LineWidth', 2);
    ylabel('Natural Frequency (Hz)')

    yyaxis right;
    plot(Qs, W_S2(:,1), 'b-', 'LineWidth', 2)
    plot(Qs, W_S2(:,2), 'r-', 'LineWidth', 2)
    plot(Qs, W_S2(:,3), 'm-', 'LineWidth', 2)

    ylabel('Sobol'' Index')
    ylim([-0.01 1.01])

    yyaxis left
    legend([aa bb cc(1)], 'Simulations', 'PCE Mean', 'PCE $\mu\pm\sigma$', 'Location', 'northeast')

    xlim(10.^[Ast Aen])
    subplot(2,1,2)
    plot(Qs, Zcofs(:,1)*100, 'k-', 'LineWidth', 2)
    plot(Qs, (Zcofs(:,1)+kron([1 -1], sqrt(Zvar)))*100, 'k-.', 'LineWidth', 2)
    ylabel('Damping Factor (\%)')
    xlabel('Modal Amplitude')
    xlim(10.^[Ast Aen])

    yyaxis right;
    aa = plot(Qs, Z_S2(:,1), 'b-', 'LineWidth', 2);
    bb = plot(Qs, Z_S2(:,2), 'r-', 'LineWidth', 2);
    cc = plot(Qs, Z_S2(:,3), 'm-', 'LineWidth', 2);
    legend([aa bb cc], ...
           sprintf('$SU_{%s,%s}$',pttls{uqparsi(Js(1,1))},pttls{uqparsi(Js(1,2))}), ...
           sprintf('$SU_{%s,%s}$',pttls{uqparsi(Js(2,1))},pttls{uqparsi(Js(2,2))}), ...
           sprintf('$SU_{%s,%s}$',pttls{uqparsi(Js(3,1))},pttls{uqparsi(Js(3,2))}), ...
           'Location', 'northeast')
    ylabel('Sobol'' Index')
    ylim([-0.01 1.01])

    if savfig
        figdir = ['./FIGS/' symnam '_Nq=' num2str(Nqpce)];
        mkdir(figdir)
        export_fig([figdir '/BB_SU2.eps'], '-depsc');
    end
end
% * Plotting PDFs
sc = 2.5;
uns = {'kN/m', 'kN/m', '', ''};

i = 2;
for i=1:Nnq
    figure(i*10)
    clf()
    set(gcf, 'Color', 'white')
    
    if min(XIs(:,i))>0
        Xs = linspace(min(XIs(:,i))/sc, max(XIs(:,i))*sc, 100);
    else
        Xs = linspace(min(XIs(:,i))*sc, max(XIs(:,i))*sc, 100);
    end
    plot(germgens(Xs,uqparsi(i)), pdffuns{uqparsi(i)}(Xs), 'LineWidth', 2); hold on
    stem(germgens(XIs(:,i),uqparsi(i)), pdffuns{uqparsi(i)}(XIs(:,i)), 'MarkerFaceColor', 'w', 'LineWidth', 1.5)
    % xlabel(sprintf('Germ $x_{%s}$', pttls{uqparsi(i)}))
    xlabel(sprintf('Factor Value $%s$ %s', pttls{uqparsi(i)}, uns{uqparsi(i)}))
    ylabel('Probability Density Function')

    if savfig
        export_fig([figdir '/PDF_' num2str(i) '.eps'], '-depsc');
    end
end

% *
fac_gIs = cell2mat(arrayfun(@(i) germgens(XIgs(:,i),uqparsi(i)), 1:Nnq, 'UniformOutput', false));
if Nnq==2
    figure(100)
    set(gcf, 'Color', 'white')
    clf()
    plot(fac_gIs(:,1), fac_gIs(:,2), 'ko', 'MarkerFaceColor', 'k');
    xlabel(sprintf('Factor Value $%s$ %s', pttls{uqparsi(1)}, uns{uqparsi(1)}))
    ylabel(sprintf('Factor Value $%s$ %s', pttls{uqparsi(2)}, uns{uqparsi(2)}))
    grid on

    if savfig
        export_fig([figdir '/PGRID.eps'], '-depsc');
    end
end

% *
load([savdir '/sim-data.mat'], 'erris');

colos = DISTINGUISHABLE_COLORS(Nqpce^(Nnq-1));
scs = [1e-3 1e-3 1.0 1.0];
uns = {'kN/m', 'kN/m', '', ''};

figure(100)
clf()
set(gcf, 'Color', 'white')
aa = gobjects(Nqpce^(Nnq-1),1);
for isim=find(~erris)'
    load(sprintf('%s/sim-%d.mat', savdir, isim), 'uxwL', 'xis', 'wis');
    UxwL = uxwL;
        
    subplot(2,1,1)
    aa(floor((isim-1)/Nqpce)+1) = semilogx(10.^UxwL(end,:), UxwL(end-1,:)/2/pi, '-', 'Color', colos(floor((isim-1)/Nqpce)+1,:)); hold on
    switch Nnq
      case 2
        legend(aa(floor((isim-1)/Nqpce)+1), ...
               sprintf('$%s$ = %f %s', pttls{uqparsi(1)}, ...
                       fac_gIs(isim,1)*scs(uqparsi(1)), uns{uqparsi(1)}))
      case 3
        legend(aa(floor((isim-1)/Nqpce)+1), ...
               sprintf('$(%s, %s)$ = (%f %s, %f %s) ', ...
                       pttls{uqparsi(1)}, pttls{uqparsi(2)}, ...
                       fac_gIs(isim,1)*scs(uqparsi(1)), uns{uqparsi(1)}, ...
                       fac_gIs(isim,2)*scs(uqparsi(2)), uns{uqparsi(2)}))
      otherwise
        warning('Nnq>3')
    end        
    
    subplot(2,1,2)
    semilogx(10.^UxwL(end,:), UxwL(end-2,:)./(2*UxwL(end-1,:))*100, '-', 'Color', colos(floor((isim-1)/Nqpce)+1,:)); hold on
    drawnow
end
subplot(2,1,1)
xlim(10.^[Ast Aen])
ll=legend(aa, 'Location', 'northeast');
ylabel('Natural Frequency (Hz)')
subplot(2,1,2)
xlim(10.^[Ast Aen])
ylabel('Damping Factor (\%)')
xlabel('Modal Amplitude')

if savfig
    if Nqpce^(Nnq-1)>6
        set(ll, 'Visible', 'off');
    end
    export_fig([figdir '/BB_data.eps'], '-depsc');
end
