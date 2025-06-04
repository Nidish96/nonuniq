%% * Preamble
% clc
clear all
addpath('./ROUTINES')
addpath('./ROUTINES/HARMONIC/')
addpath('./ROUTINES/SOLVERS/')
addpath('./ROUTINES/QUADRATURE/')
addpath('./ROUTINES/export_fig/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

analyze = false;
savfigs = true;

%% * Parameters (Modified, from Erhan)
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

%% * Setup Friction Model Parameters
fpars = [deg2rad(40);3e5;3e5;0.9;0]; % [alpha; kt;kn;mu;f0]

nsig = 1;
kvals = 3e5*(1+[-nsig 0 nsig]*.3);
parv = 'kn';

mpars = struct('M', M, 'C', C, 'K', K, ...
               'Lrels', @(al) cat(3, Lrel1(al), Lrel2(al)));

%% * Setup HB
h = 0:5;
Nhc = sum((h==0)+2*(h~=0));
[inds0, indsh, zinds, rinds, iinds] = HINDS(6, h(:));
Nt = 256;

Fstat = zeros(6*Nhc,1);
Fstat(zinds(1:6)) = F0;

mi = 3; % Mode of interest

UxwL = cell(size(kvals));
Lscs = zeros(1, length(kvals));
if analyze
    for ik=1:length(kvals)
        switch parv
          case 'kt'
            fpars(2) = kvals(ik);
          case 'kn'
            fpars(3) = kvals(ik);            
        end

        % Linear Eigenanalysis
        K0 = mpars.K + Lrel1(fpars(1))'*diag([fpars(2) fpars(3)])*Lrel1(fpars(1)) ...
             + Lrel2(fpars(1))'*diag([fpars(2) fpars(3)])*Lrel2(fpars(1));
        mpars.dK = K0-mpars.K;
        [V, Wst] = eig(K0, M);
        [Wst, si] = sort(sqrt(diag(Wst)));
        V = V(:,si);
        V = V./sqrt(diag(V'*M*V)');

        % Initial Guess
        Uxw0 = zeros(6*Nhc+2,1);
        Uxw0(zinds) = K0\F0;
        % Uxw0([rinds(1:6) iinds(1:6)]) = [V(:, mi); V(:, mi)]/sqrt(2);
        Uxw0(iinds(1:6)) = V(:,mi);
        Uxw0(end-1) = V(:,mi)'*C*V(:,mi);
        Uxw0(end) = Wst(mi);

        % Obtain length scaling
        utunstat = [Lrel1(fpars(1))*(K0\F0) Lrel2(fpars(1))*(K0\F0)];
        utunms = [Lrel1(fpars(1))*V(:,mi) Lrel2(fpars(1))*V(:,mi)];
        
        fn0 = fpars(3)*utunstat(2,:);
        fnm = fpars(3)*utunms(2,:);

        ft0 = fpars(2)*utunstat(1,:);
        ftm = fpars(2)*utunms(1,:);
        ftm = ftm.*sign(ft0);
        ft0 = ft0.*sign(ft0);

        slims = fpars(4)*fpars(3)*utunstat(2,:);  % mu*kn*un0 -> fslip

        aslips = slims./(fpars(2)*utunms(1,:));

        aslips = slims./(fpars(2)*abs(utunms(1,:))-...
		         abs(fpars(4)*fpars(3)*utunms(2,:)));

        aslips = (fpars(4)*fn0)./(abs(ftm)+fpars(4)*abs(fnm));

        aslips = (fpars(4)*fn0)./(abs(ftm)-fpars(4)*abs(fnm));

        length_scale = min(abs(aslips));
        Lscs(ik) = length_scale;

        % Construct matrix for implicitly enforcing phase constraint
        Lb = speye(Nhc*6+2);	Lb(:, rinds(5)) = [];

        Ast = -4;  % -4
        Aen = -1.5;

        Ast = log10(1e-1*length_scale);
        Aen = log10(1e2*length_scale);
        ds = 0.05; % 0.1
                   % ds = 0.02; % [0.25, 0.1, 0.2 0.5 0.05]
        DSS = [0.25 0.1 0.2 0.5 0.05 0.02 0.06 0.025];

        Sopt = struct('stepmax', 5000, 'jac', 'xl', 'dsmin', 0.001, 'dynamicDscale', 1);
        Sopt.Dscale = [ones(6*Nhc-1,1)*1e-4; (V(:,mi)'*C*V(:,mi)); Wst(mi); 1];

        flg = false;
        
        % Forward continuation (this is preferred always)
        for i=1:length(DSS)
            ds = DSS(i);
            UxwL{ik} = solve_and_continue(Lb'*Uxw0, ...
                                          @(Uxwl) EPMCRESFUN(Uxwl, Fstat, Lb, h, Nt, fpars, mpars), ...
                                          Ast, Aen, ds, Sopt);
            if UxwL{ik}(end)>Aen
                flg = true;
                break
            end
        end
        if ~flg % Try in reverse
            for i=1:length(DSS)
                ds = DSS(i);    
                UxwL{ik} = solve_and_continue(Lb'*Uxw0, ...
                                              @(Uxwl) EPMCRESFUN(Uxwl, Fstat, Lb, h, ...
                                                                 Nt, fpars, mpars), ...
                                              Aen, Ast, ds, Sopt);
                if UxwL{ik}(end)<Ast
                    UxwL{ik} = UxwL{ik}(:, end:-1:1);
                    flg = true;
                    break
                end    
            end
        end

        if ~flg
            warning('Incomplete Continuation')
        end
        UxwL{ik} = [Lb*UxwL{ik}(1:end-1,:); UxwL{ik}(end,:)];
    end

    save(sprintf('./DATS/C3_UxwL_%s.mat', parv), 'UxwL', 'Lscs', 'kvals', 'fpars');
else
    load(sprintf('./DATS/C3_UxwL_%s.mat', parv), 'UxwL', 'Lscs', 'kvals', 'fpars');
end

%% * Plot Backbones
colos = DISTINGUISHABLE_COLORS(length(kvals));
psymb = '-';
fsz = 16;
xls = [1e0 1e1];
% xls = (10.^[Ast Aen])/length_scale;

Asts = log10(1e-1*Lscs);
Aens = log10(1e2*Lscs);
xls = 10.^[floor(min(Asts)) ceil(max(Aens))];
xls(2) = 1e-1;

figure(1);
poss=get(gcf, 'Position');
set(gcf, 'Color', 'white', 'Position', [poss(1:2) 920 420])
clf()
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact')

ax1=nexttile;
aa = gobjects(size(kvals));
for ik=1:length(kvals)
    aa(ik) = semilogx(10.^UxwL{ik}(end,:), ...
                      UxwL{ik}(end-1,:)/2/pi, psymb, 'LineWidth', 2, ...
                      'Color', colos(ik,:)); hold on
    switch parv
      case 'kt'
        legend(aa(ik), sprintf('$k_t=%.2f\\,\\mathrm{MN\\,m^{-1}}$', kvals(ik)*1e-6))
      case 'kn'
        legend(aa(ik), sprintf('$k_n=%.2f\\,\\mathrm{MN\\,m^{-1}}$', kvals(ik)*1e-6))
    end
end
legend(aa, 'Location', 'southwest');
grid on; grid minor; grid minor 
% plot(length_scale*[1 1], ylim, '--')
ylabel('Natural Frequency (Hz)')
xlabel('Modal Amplitude')
xlim(xls);
set(gca, 'FontSize', fsz, 'XTick', 10.^(-5:2:-1))

ax2=nexttile;
for ik=1:length(kvals)
    semilogx(10.^UxwL{ik}(end,:), ...
             UxwL{ik}(end-2,:)./(2*UxwL{ik}(end-1,:))*100, psymb, 'LineWidth', 2, ...
             'Color', colos(ik,:));
    hold on
end
grid on; grid minor; grid minor 
% plot(length_scale*[1 1], ylim, '--')
ylabel('Damping Factor (\%)')
xlabel('Modal Amplitude')
xlim(xls);
set(gca, 'FontSize', fsz, 'XTick', 10.^(-5:2:-1))

if savfigs
    export_fig(sprintf('./FIGS/C3_DEMOBB_0_%s.eps',parv), '-depsc');
end

%% ** Plot Scaled
figure(2);
poss=get(gcf, 'Position');
set(gcf, 'Color', 'white', 'Position', [poss(1:2) 920 420])
clf()
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact')

ax1=nexttile;
aa = gobjects(size(kvals));
for ik=1:length(kvals)
    aa(ik) = semilogx(10.^UxwL{ik}(end,:)/Lscs(ik), ...
                      UxwL{ik}(end-1,:)/2/pi, psymb, 'LineWidth', 2, ...
                      'Color', colos(ik,:)); hold on
    switch parv
      case 'kt'
        legend(aa(ik), sprintf('$k_t=%.2f\\,\\mathrm{MN\\,m^{-1}}$', kvals(ik)*1e-6))
      case 'kn'
        legend(aa(ik), sprintf('$k_n=%.2f\\,\\mathrm{MN\\,m^{-1}}$', kvals(ik)*1e-6))
    end
end
legend(aa, 'Location', 'southwest');
grid on; grid minor; grid minor 
% plot(length_scale*[1 1], ylim, '--')
ylabel('Natural Frequency (Hz)')
xlabel('Scaled Modal Amplitude')
set(gca, 'FontSize', fsz)
xlim([10^(-1.5) 1e2])

ax2=nexttile;
for ik=1:length(kvals)
    semilogx(10.^UxwL{ik}(end,:)/Lscs(ik), ...
             UxwL{ik}(end-2,:)./(2*UxwL{ik}(end-1,:))*100, psymb, 'LineWidth', 2, ...
             'Color', colos(ik,:));
    hold on
end
grid on; grid minor; grid minor 
% plot(length_scale*[1 1], ylim, '--')
ylabel('Damping Factor (\%)')
xlabel('Scaled Modal Amplitude')
set(gca, 'FontSize', fsz)
xlim([10^(-1.5) 1e2])

if savfigs
    export_fig(sprintf('./FIGS/C3_DEMOBB_2_%s.eps',parv), '-depsc');
end

%% * Plot for paper
colos = DISTINGUISHABLE_COLORS(length(kvals));
psymb = '-';
fsz = 20;
xls = [1e0 1e1];
% xls = (10.^[Ast Aen])/length_scale;

Asts = log10(1e-1*Lscs);
Aens = log10(1e2*Lscs);
xls = 10.^[floor(min(Asts)) ceil(max(Aens))];
xls(2) = 1e-1;

figure(10);
poss=get(gcf, 'Position');
set(gcf, 'Color', 'white', 'Position', poss)
clf()
tiledlayout(1,1, 'TileSpacing', 'compact', 'Padding', 'compact')

ax1=nexttile;
aa = gobjects(size(kvals));
for ik=1:length(kvals)
    aa(ik) = semilogx(10.^UxwL{ik}(end,:)/Lscs(ik), ...
                      UxwL{ik}(end-1,:)/2/pi, psymb, 'LineWidth', 2, ...
                      'Color', colos(ik,:)); hold on
    switch parv
      case 'kt'
        legend(aa(ik), sprintf('$k_t=%.2f\\,\\mathrm{MN\\,m^{-1}}$', kvals(ik)*1e-6))
      case 'kn'
        legend(aa(ik), sprintf('$k_n=%.2f\\,\\mathrm{MN\\,m^{-1}}$', kvals(ik)*1e-6))
    end
end
legend(aa, 'Location', 'southwest');

grid on; grid minor; grid minor 
% plot(length_scale*[1 1], ylim, '--')
ylabel('Frequency [Hz]')
xlabel('Scaled Modal Amplitude $q/\ell$')
% xlim(xls);
set(gca, 'FontSize', fsz)

figure(11);
poss=get(gcf, 'Position');
set(gcf, 'Color', 'white', 'Position', poss)
clf()
tiledlayout(1,1, 'TileSpacing', 'compact', 'Padding', 'compact')

ax2=nexttile;
for ik=1:length(kvals)
    semilogx(10.^UxwL{ik}(end,:)/Lscs(ik), ...
             UxwL{ik}(end-2,:)./(2*UxwL{ik}(end-1,:))*100, psymb, 'LineWidth', 2, ...
             'Color', colos(ik,:));
    hold on
end
grid on; grid minor; grid minor 
% plot(length_scale*[1 1], ylim, '--')
ylabel('Damping Factor [\%]')
xlabel('Scaled Modal Amplitude $q/\ell$')
% xlim(xls);
set(gca, 'FontSize', fsz)

ax = axes('Position', [0.26 0.58 0.18 0.3]);
bb = bar(kvals*1e-6, Lscs*1e3, 'FaceColor', 'flat');
bb.CData = colos;
bb.Labels = arrayfun(@(l) sprintf('%.2f', l*1e3), Lscs, 'UniformOutput', false);
switch parv
  case 'kt'
    xlabel('$k_t$ [$\mathrm{MN\,m^{-1}}$]')
  case 'kn'
    xlabel('$k_n$ [$\mathrm{MN\,m^{-1}}$]')
end
ylabel('$\ell$ [$\mathrm{kg^{1/2}\,mm}$]')
set(gca, 'FontSize', fix(fsz*0.75))
grid on;grid minor;grid minor;box on;

if savfigs
    figure(10)
    export_fig(sprintf('./FIGS/C3_WBB_%s.pdf', parv), '-dpdf')

    figure(11)
    export_fig(sprintf('./FIGS/C3_ZBB_%s.pdf', parv), '-dpdf')    
end
