%% * Preamble
% clc
clear all
addpath('./ROUTINES')
addpath('./ROUTINES/HARMONIC/')
addpath('./ROUTINES/SOLVERS/')
addpath('./ROUTINES/QUADRATURE/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

analyze = false;
savfigs = false;

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
fpars = [deg2rad(40);3e5;3e5;0.9;-1]; % [alpha; kt;kn;mu;f0]

mvals = -1:.5:1;

mpars = struct('M', M, 'C', C, 'K', K, ...
               'Lrels', @(al) cat(3, Lrel1(al), Lrel2(al)));

%% * Linear Eigen Analysis (stick)
K0 = mpars.K + Lrel1(fpars(1))'*diag([fpars(2) fpars(3)])*Lrel1(fpars(1)) ...
     + Lrel2(fpars(1))'*diag([fpars(2) fpars(3)])*Lrel2(fpars(1));
mpars.dK = K0-mpars.K;
[V, Wst] = eig(K0, M);
[Wst, si] = sort(sqrt(diag(Wst)));
V = V(:,si);
V = V./sqrt(diag(V'*M*V)');

%% * Setup HB
h = 0:5;
Nhc = sum((h==0)+2*(h~=0));
[inds0, indsh, zinds, rinds, iinds] = HINDS(6, h(:));
Nt = 256;

Fstat = zeros(6*Nhc,1);
Fstat(zinds(1:6)) = F0;

mi = 3; % Mode of interest

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

UxwL = cell(size(mvals));
if analyze
    for im=1:length(mvals)
        fpars(5) = mvals(im);
        
        % Forward continuation (this is preferred always)
        for i=1:length(DSS)
                 ds = DSS(i);
        UxwL{im} = solve_and_continue(Lb'*Uxw0, ...
                                      @(Uxwl) EPMCRESFUN(Uxwl, Fstat, Lb, h, Nt, fpars, mpars), ...
                                      Ast, Aen, ds, Sopt);
        if UxwL{im}(end)>Aen
            flg = true;
            break
        end
    end
        if ~flg % Try in reverse
                 for i=1:length(DSS)
            ds = DSS(i);    
            UxwL{im} = solve_and_continue(Lb'*Uxw0, ...
                                          @(Uxwl) EPMCRESFUN(Uxwl, Fstat, Lb, h, ...
                                                             Nt, fpars, mpars), ...
                                          Aen, Ast, ds, Sopt);
            if UxwL{im}(end)<Ast
                UxwL{im} = UxwL{im}(:, end:-1:1);
                flg = true;
                break
            end    
        end
    end

        if ~flg
            warning('Incomplete Continuation')
        end
        UxwL{im} = [Lb*UxwL{im}(1:end-1,:); UxwL{im}(end,:)];
    end

    save('./DATS/C1_UxwL.mat', 'UxwL', 'mvals', 'fpars');
else
    load('./DATS/C1_UxwL.mat', 'UxwL', 'mvals', 'fpars');
end

%% * Plot Backbones
colos = DISTINGUISHABLE_COLORS(length(mvals));
psymb = '-';
fsz = 16;
xls = [1e0 1e1];
% xls = (10.^[Ast Aen])/length_scale;

figure(1);
poss=get(gcf, 'Position');
set(gcf, 'Color', 'white', 'Position', [poss(1:2) 920 420])
clf()
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact')

ax1=nexttile;
aa = gobjects(size(mvals));
for im=1:length(mvals)
    aa(im) = semilogx(10.^UxwL{im}(end,:)/length_scale, ...
                      UxwL{im}(end-1,:)/2/pi, psymb, 'LineWidth', 2, ...
                      'Color', colos(im,:)); hold on
    legend(aa(im), sprintf('$m=%.1f$', mvals(im)))
end
legend(aa, 'Location', 'southwest');
grid on; grid minor; grid minor 
% plot(length_scale*[1 1], ylim, '--')
ylabel('Natural Frequency (Hz)')
xlabel('Scaled Modal Amplitude')
set(gca, 'FontSize', fsz)

ax2=nexttile;
for im=1:length(mvals)
    semilogx(10.^UxwL{im}(end,:)/length_scale, ...
             UxwL{im}(end-2,:)./(2*UxwL{im}(end-1,:))*100, psymb, 'LineWidth', 2, ...
             'Color', colos(im,:));
    hold on
end
grid on; grid minor; grid minor 
% plot(length_scale*[1 1], ylim, '--')
ylabel('Damping Factor (\%)')
xlabel('Scaled Modal Amplitude')
set(gca, 'FontSize', fsz)

if savfigs
    export_fig('./FIGS/C1_DEMOBB_0.eps', '-depsc');

    for ax=[ax1 ax2]
        axes(ax);
        xlim(xls);
    end
    export_fig('./FIGS/C1_DEMOBB_1.eps', '-depsc');
end

%% ** Plot for paper
fsz = 20;

figure(10);
poss=get(gcf, 'Position');
set(gcf, 'Color', 'white', 'Position', poss)
clf()
tiledlayout(1,1, 'TileSpacing', 'compact', 'Padding', 'compact')

ax1=nexttile;
aa = gobjects(size(mvals));
for im=1:length(mvals)
    aa(im) = semilogx(10.^UxwL{im}(end,:)/length_scale, ...
                      UxwL{im}(end-1,:)/2/pi, psymb, 'LineWidth', 2, ...
                      'Color', colos(im,:)); hold on
    legend(aa(im), sprintf('$m=%.1f$', mvals(im)))
end
legend(aa, 'Location', 'southwest');
grid on; grid minor; grid minor 
% plot(length_scale*[1 1], ylim, '--')
ylabel('Frequency [Hz]')
xlabel('Scaled Modal Amplitude $q/\ell$')
set(gca, 'FontSize', fsz)

figure(11);
poss=get(gcf, 'Position');
set(gcf, 'Color', 'white', 'Position', poss)
clf()
tiledlayout(1,1, 'TileSpacing', 'compact', 'Padding', 'compact')

ax2=nexttile;
for im=1:length(mvals)
    semilogx(10.^UxwL{im}(end,:)/length_scale, ...
             UxwL{im}(end-2,:)./(2*UxwL{im}(end-1,:))*100, psymb, 'LineWidth', 2, ...
             'Color', colos(im,:));
    hold on
end
grid on; grid minor; grid minor 
% plot(length_scale*[1 1], ylim, '--')
ylabel('Damping Factor [\%]')
xlabel('Scaled Modal Amplitude $q/\ell$')
set(gca, 'FontSize', fsz)

if savfigs
    figure(10)
    export_fig(sprintf('./FIGS/C1_WBB.pdf'), '-dpdf')

    figure(11)
    export_fig(sprintf('./FIGS/C1_ZBB.pdf'), '-dpdf')    
end


%% * Depict Mode Shape
ell1 = 1;
ell2 = 0.8;
angL = (pi-fpars(1))/2;

dell = (ell1-ell2);

obL = [0 0;
       ell2 0;
       ell1 dell*tan(angL);
       ell1 ell1;
       0 ell1];
wd = [-dell 0;
      0 dell*tan(angL);
      dell 0];

ob1 = [obL(:,1)-1, obL(:,2)];
ob2 = [1-obL(:,1), obL(:,2)];

ob1 = (ob1-mean(ob1))*0.9+mean(ob1);
ob2 = (ob2-mean(ob2))*0.9+mean(ob2);
wd  = (wd-mean(wd))*0.9+mean(wd);


V0 = UxwL{1}(rinds(1:6),1)+1j*UxwL{1}(iinds(1:6),1); V0 = V0*exp(-1j*angle(V0(1)));
V8 = UxwL{1}(rinds(1:6),end)+1j*UxwL{1}(iinds(1:6),end); V8 = V8*exp(-1j*angle(V8(1)));

Vp = V(:,mi);

% Vp = real(V8);

% mi = 3;
sc = 0.25;

figure(2);
set(gcf, 'Color', 'white')
clf()
% plot(ob1([1:end 1],1), ob1([1:end 1],2), 'b-.', 'LineWidth', 2); hold on
% plot(ob2([1:end 1],1), ob2([1:end 1],2), 'r-.', 'LineWidth', 2); hold on
% plot(wd([1:end 1],1), wd([1:end 1],2), 'g-.', 'LineWidth', 2); hold on

fill(ob1(:,1), ob1(:,2), 'w', 'LineWidth', 2, 'FaceAlpha', 0.6, 'EdgeColor', 'b', 'LineStyle', '-.'); hold on
fill(ob2(:,1), ob2(:,2), 'w', 'LineWidth', 2, 'FaceAlpha', 0.6, 'EdgeColor', 'r', 'LineStyle', '-.');
fill(wd(:,1), wd(:,2), 'w', 'LineWidth', 2, 'FaceAlpha', 0.6, 'EdgeColor', 'g', 'LineStyle', '-.');

fill(ob1(:,1)+Vp(1)*sc, ob1(:,2)+Vp(2)*sc, ...
     'b', 'LineWidth', 2, 'FaceAlpha', 0.6); hold on
fill(ob2(:,1)+Vp(3)*sc, ob2(:,2)+Vp(4)*sc, ...
     'r', 'LineWidth', 2, 'FaceAlpha', 0.6)
fill(wd(:,1)+Vp(5)*sc, wd(:,2)+Vp(6)*sc, ...
     'g', 'LineWidth', 2, 'FaceAlpha', 0.6);

fill(ob1(:,1)-Vp(1)*sc, ob1(:,2)-Vp(2)*sc, ...
     'b', 'LineWidth', 2, 'FaceAlpha', 0.2, ...
     'EdgeAlpha', 0.2); hold on
fill(ob2(:,1)-Vp(3)*sc, ob2(:,2)-Vp(4)*sc, ...
     'r', 'LineWidth', 2, 'FaceAlpha', 0.2, ...
     'EdgeAlpha', 0.2)
fill(wd(:,1)-Vp(5)*sc, wd(:,2)-Vp(6)*sc, ...
     'g', 'LineWidth', 2, 'FaceAlpha', 0.2, ...
     'EdgeAlpha', 0.2);

axis equal
axis off
grid on
xlim([-0.8 1]*1.25)
ylim([-0.25 1]*1.5)

if savfigs
    export_fig(sprintf('./FIGS/C1_DEMOMSHAPE_%d.eps', mi), '-depsc')
end
