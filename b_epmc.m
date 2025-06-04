clc
clear all
addpath('./ROUTINES/HARMONIC/')
addpath('./ROUTINES/SOLVERS/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13);

%% Parameters
pars = struct('m', 1, 'mD', 0.1, 'c', 20, ...
    'k', 3e5, 'kxy', 3e5, 'k12', 7e5, 'kD', 3e5, ...
    'kn', 3e5, 'kt', 3e5, 'mu', 0.5);

M = diag([repmat(pars.m, 1,4) repmat(pars.mD, 1,2)]);
K = diag([repmat(pars.k, 1,4) zeros(1,2)]);
K([1 3], [1 3]) = K([1 3], [1 3]) + pars.k12*[1 -1;-1 1];
K([2 4], [2 4]) = K([2 4], [2 4]) + pars.k12*[1 -1;-1 1];
K(1:2,1:2) = K(1:2,1:2) + pars.kxy*[1 -1;-1 1];
K(3:4,3:4) = K(3:4,3:4) + pars.kxy*[1 -1;-1 1];
K(5:6,5:6) = K(5:6,5:6) + pars.kD*[1 -1;-1 1];
C = diag([repmat(pars.c, 1,4) zeros(1,2)]);
F0 = [0; 0; 0; 0; 0; 120];
Fv = [0; -20; 5; 0; 0; 0];

al = deg2rad(30);  % alpha (angle value)
Lrel1 = @(al) kron([-1 0 1], [cos(al) sin(al);-sin(al) cos(al)]);
Lrel2 = @(al) kron([0 -1 1], [cos(al) -sin(al);sin(al) cos(al)]);

mpars = struct('M', M, 'C', C, 'K', K, ...
    'Lrels', @(al) cat(3, Lrel1(al), Lrel2(al)));
fpars = [deg2rad(30);pars.kt;pars.kn;pars.mu;0]; % [alpha; kt;kn;mu;f0]

K0 = mpars.K + pars.kt*(Lrel1(fpars(1))'*Lrel1(fpars(1))) +  pars.kt*(Lrel2(fpars(1))'*Lrel2(fpars(1)));
mpars.K0 = K0;

%% Linear Eigen Analysis (stick)
[V, Wst] = eig(K0, M);
[Wst, si] = sort(sqrt(diag(Wst)));
V = V(:,si);

%% Check Residue
h = 0:5;
Nhc = sum((h==0)+2*(h~=0));
[inds0, indsh, zinds, rinds, iinds] = HINDS(6, h(:));
Nt = 256;

Fstat = zeros(6*Nhc,1);
Fstat(zinds(1:6)) = F0;

Fl = zeros(6*Nhc,1);
Fl(rinds(5)) = 1.0;

%%
mi = 3;

Uxw0 = zeros(6*Nhc+2,1);
Uxw0(zinds) = K0\F0;
Uxw0([rinds(1:6) iinds(1:6)]) = [V(:, mi); V(:, mi)]/sqrt(2);
Uxw0(end-1) = V(:,mi)'*C*V(:,mi);
Uxw0(end) = Wst(mi);

Ast = -5;
Aen = 0;
ds = 0.1;
fpars = [deg2rad(30);pars.kt;pars.kn;pars.mu;1]; % [alpha; kt;kn;mu;f0]

opt = struct('Display', true);
UxwS = NSOLVE(@(Uxw) EPMCRESFUN([Uxw;Ast], Fstat, Fl, h, Nt, fpars, mpars), Uxw0, opt);

% opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% UxwS = fsolve(@(Uxw) EPMCRESFUN([Uxw;-4], Fstat, Fl, h, Nt, fpars, mpars), Uxw0, opt);

% Copt = struct('Nmax', 1000);
% UxwL = CONTINUE(@(Uxwl) EPMCRESFUN(Uxwl, Fstat, Fl, h, Nt, fpars, mpars), Uxw0, Ast, Aen, ds, Copt);

Sopt = struct('stepmax', 2000, 'jac', 'xl');
UxwL = solve_and_continue(Uxw0, @(Uxwl) EPMCRESFUN(Uxwl, Fstat, Fl, h, Nt, fpars, mpars), Ast, Aen, ds, Sopt);

figure(2);
% clf()
set(gcf, 'Color', 'white')
subplot(2,1,1)
semilogx(10.^UxwL(end,:), UxwL(end-1,:)/2/pi, '.-'); hold on
ylabel('Natural Frequency (Hz)')

subplot(2,1,2)
semilogx(10.^UxwL(end,:), UxwL(end-2,:)./(2*UxwL(end-1,:))*100, '.-'); hold on
ylabel('Damping Factor (\%)')
xlabel('Modal Amplitude')
