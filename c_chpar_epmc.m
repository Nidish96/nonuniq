% clc
clear all
addpath('./ROUTINES/HARMONIC/')
addpath('./ROUTINES/SOLVERS/')
addpath('./ROUTINES/QUADRATURE/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

%% Parameters (Modified, from Erhan)
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

%% Setup Friction Model Parameters
fpars = [deg2rad(40);3e5;3e5;0.6;-1]; % [alpha; kt;kn;mu;f0]

% fpars = [deg2rad(40);3e5;3e5;1.8;-1]; % [alpha; kt;kn;mu;f0]
% fpars(2:5) = [2.2997e5 3.2226e5 0.5 0.33998];
% fpars(2:5) = [3e5 3e5 0.3 0];

% fpars(4) = 1.6;

% al = fpars(1);
mpars = struct('M', M, 'C', C, 'K', K, ...
               'Lrels', @(al) cat(3, Lrel1(al), Lrel2(al)));

%% Linear Eigen Analysis (stick)
K0 = mpars.K + Lrel1(fpars(1))'*diag([fpars(2) fpars(3)])*Lrel1(fpars(1)) ...
     + Lrel2(fpars(1))'*diag([fpars(2) fpars(3)])*Lrel2(fpars(1));
mpars.dK = K0-mpars.K;
[V, Wst] = eig(K0, M);
[Wst, si] = sort(sqrt(diag(Wst)));
V = V(:,si);
V = V./sqrt(diag(V'*M*V)');

%% Setup HB
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

length_scale = min(abs(aslips))

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
    UxwL = solve_and_continue(Lb'*Uxw0, ...
                              @(Uxwl) EPMCRESFUN(Uxwl, Fstat, Lb, h, Nt, fpars, mpars), ...
                              Ast, Aen, ds, Sopt);
    if UxwL(end)>Aen
        flg = true;
        break
    end
end
if ~flg % Try in reverse
    for i=1:length(DSS)
        ds = DSS(i);    
        UxwL = solve_and_continue(Lb'*Uxw0, ...
                                  @(Uxwl) EPMCRESFUN(Uxwl, Fstat, Lb, h, ...
                                                     Nt, fpars, mpars), ...
                                  Aen, Ast, ds, Sopt);
        if UxwL(end)<Ast
            UxwL = UxwL(:, end:-1:1);
            flg = true;
            break
        end    
    end
end

if ~flg
    warning('Incomplete Continuation')
end
UxwL = [Lb*UxwL(1:end-1,:); UxwL(end,:)];

%%
figure(1);
% clf()
set(gcf, 'Color', 'white')
tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile;
semilogx(10.^UxwL(end,:)/length_scale, UxwL(end-1,:)/2/pi, '.-'); hold on
grid on; grid minor; grid minor 
% plot(length_scale*[1 1], ylim, '--')
ylabel('Natural Freq. (Hz)')

nexttile;
semilogx(10.^UxwL(end,:)/length_scale, UxwL(end-2,:)./(2*UxwL(end-1,:))*100, '.-');
hold on
grid on; grid minor; grid minor 
% plot(length_scale*[1 1], ylim, '--')
ylabel('Damping Factor (\%)')
xlabel('Scaled Modal Amplitude')
