clc
clear all
addpath('./ROUTINES/HARMONIC/')
addpath('./ROUTINES/SOLVERS/')

%% Parameters
pars = struct('m', 1, 'mD', 0.1, 'c', 20, ...
    'k', 3e5, 'kxy', 3e5, 'k12', 7e5, 'kD', 3e5, ...
    'kn', 3e5, 'kt', 3e5, 'mu', 0.5);
% Hellö
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

%% Check Residue
h = 0:5;
Nhc = sum((h==0)+2*(h~=0));
[inds0, indsh, zinds, rinds, iinds] = HINDS(6, h(:));
Nt = 256;

mpars = struct('M', M, 'C', C, 'K', K, ...
    'Lrels', @(al) cat(3, Lrel1(al), Lrel2(al)));
fpars = [deg2rad(30);pars.kt;pars.kn;pars.mu;0]; % [alpha; kt;kn;mu;f0]

K0 = mpars.K + pars.kt*(Lrel1(fpars(1))'*Lrel1(fpars(1))) +  pars.kt*(Lrel2(fpars(1))'*Lrel2(fpars(1)));

Fl = zeros(6*Nhc,1);
Fl(zinds(1:6)) = F0;
Fl(iinds(1:6)) = Fv;

Wst = 2*pi*160;
K0 = mpars.K + pars.kt*(Lrel1(fpars(1))'*Lrel1(fpars(1))) +  pars.kt*(Lrel2(fpars(1))'*Lrel2(fpars(1)));
U0 = HARMONICSTIFFNESS(mpars.M, mpars.C, K0, Wst, h)\Fl;

opt = struct('Display', true);
Us = NSOLVE(@(U) HBRESFUN([U;Wst], Fl, h, Nt, fpars, mpars), U0, opt);

% opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% Us = fsolve(@(U) HBRESFUN([U;Wst], Fl, h, Nt, fpars, mpars), U0, opt);

%%
Copt = struct('Nmax', 500);

Wst = 145*2*pi;
Wen = 235*2*pi;
ds = 0.5;

fpars = [deg2rad(30);pars.kt;pars.kn;pars.mu;0]; % [alpha; kt;kn;mu;m]
UwC = CONTINUE(@(Uw) HBRESFUN(Uw, Fl, h, Nt, fpars, mpars), Us, Wst, Wen, ds, Copt);

figure(1)
% clf()
hold on
plot(UwC(end,:)/2/pi, abs(UwC(rinds(4),:)+1j*UwC(iinds(4),:)), '.-')
xlabel('Frequency (Hz)')
ylabel('Response Amplitude (m)')
ylim([0 3.5e-4])
