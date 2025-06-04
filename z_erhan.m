clear all;
clc;
format long;
% % close all
global D Fext kt kn mu trials counter T1GtoL T1LtoG T2GtoL T2LtoG f n u

%% Damper Angles
Beta=120;

beta=.5*Beta/180*pi;
alfa=pi/2-beta;

%% Rotation Matrices
T1GtoL = [cos(alfa) sin(alfa);-sin(alfa) cos(alfa)];
T1LtoG = T1GtoL';

T2GtoL = [cos(alfa) -sin(alfa);sin(alfa) cos(alfa)];
T2LtoG = T2GtoL';

%% Contact Element Propoerties
kn = 3e5;
kt = kn;
mu = 0.5;

%% Linear System Properties
kX = 3e5;kY=3e5;
kcross=3e5;kcross2=7e5;
m = 1;md=0.1;
c = 20;

%% Forcing
Fx1s = 5;Fx1c=0;
Fy1s = 0;Fy1c=0;
Fx2s = 0;Fx2c=0;
Fy2s = -20;Fy2c=0;

%% Frequency Range
om =235*2*pi:-1:145*2*pi;

%% Free and Fully Stuck Linear Systems

%with damper
M = [m 0 0 0 0 0;0 m 0 0 0 0;0 0 md 0 0 0;0 0 0 md 0 0;0 0 0 0 m 0;0 0 0 0 0 m];

K = [(kX+kcross+kcross2) -kcross 0 0 -kcross2 0;-kcross (kY+kcross+kcross2) 0 0 0 -kcross2;0 0 kcross -kcross 0 0;0 0 -kcross kcross 0 0;...
    -kcross2 0 0 0 (kX+kcross+kcross2) -kcross; 0 -kcross2 0 0 -kcross (kY+kcross+kcross2)];

C = [c 0 0 0 0 0;0 c 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0; 0 0 0 0 c 0;0 0 0 0 0 c];


%without damper
Mfree = [m 0 0 0;0 m 0 0;0 0 m 0;0 0 0 m];

Kfree = [(kX+kcross+kcross2) -kcross -kcross2 0;-kcross (kY+kcross+kcross2) 0 -kcross2;...
    -kcross2 0 (kX+kcross+kcross2) -kcross; 0 -kcross2 -kcross (kY+kcross+kcross2)];

Cfree = [c 0 0 0;0 c 0 0;0 0 c 0;0 0 0 c];

%rotation of contact elements from local to global

ks1L = [kt;kn];
ks1G = T1LtoG * ks1L;
Ks1G = [ks1G(1) 0 -ks1G(1) 0;0 ks1G(2) 0 -ks1G(2);-ks1G(1) 0 ks1G(1) 0; 0 -ks1G(2) 0 ks1G(2)];

ks2L = [kt;kn];
ks2G = T2LtoG * ks2L;
Ks2G = [ks2G(1) 0 -ks2G(1) 0;0 ks2G(2) 0 -ks2G(2);-ks2G(1) 0 ks2G(1) 0; 0 -ks2G(2) 0 ks2G(2)];

Ks1G_11=Ks1G(1:2,1:2);Ks1G_12 = Ks1G(1:2,3:4);Ks1G_21=Ks1G(3:4,1:2); Ks1G_22=Ks1G(3:4,3:4);
Ks2G_11=Ks2G(1:2,1:2);Ks2G_12 = Ks2G(1:2,3:4);Ks2G_21=Ks2G(3:4,1:2); Ks2G_22=Ks2G(3:4,3:4);
KstuckRot = [Ks1G_11,Ks1G_12,zeros(2,2);Ks1G_21,Ks1G_22+Ks2G_11,Ks2G_12;zeros(2,2),Ks2G_21,Ks2G_22];

Ks = K + KstuckRot;

[Phi_lin,Wn_lin]=eig(K,M);
[Phi_stuck,Wn_stuck]=eig(Ks,M);
[Phi_linfree,Wn_linfree]=eig(Kfree,Mfree);
Fext1 = 1i*[ Fx1s;Fy1s;0;0;Fx2s;Fy2s];
Fextfree = 1i*[ Fx1s;Fy1s;Fx2s;Fy2s];
% Linear response.
Ylin = zeros(6,length(om));
Ylinfree = zeros(4,length(om));
Ystick = zeros(6,length(om));
for nw = 1:length(om)
    
    % Frequency.
    w = om(nw);
    
    % Free response.
    Ylin(:,nw) = (K- w^2*M + 1i*w*C)\ Fext1;
    Ylinfree(:,nw) = (Kfree- w^2*Mfree + 1i*w*Cfree)\ Fextfree;
    % Stuck response.
    Ystick(:,nw) = (Ks- w^2*M + 1i*w*C)\ Fext1;
    
end

%  Plots.
figure(10);
hold on
% plot(om/2/pi,abs(Ylin(1,:)),'g','linewidth',2);
% plot(om/2/pi,abs(Ylin(2,:)),'m','linewidth',2);
% plot(om/2/pi,abs(Ylin(5,:)),'r','linewidth',2);
plot(om/2/pi,abs(Ylin(6,:)),'B','linewidth',2);
% hold on;
% plot(om/2/pi,abs(Ystick(1,:)),'g-.','linewidth',2);
% plot(om/2/pi,abs(Ystick(2,:)),'c-.','linewidth',2);
% plot(om/2/pi,abs(Ystick(5,:)),'r--','linewidth',2);
plot(om/2/pi,abs(Ystick(6,:)),'k--','linewidth',2);
set(gca, 'YScale', 'lin')
xlim([om(end)/2/pi om(1)/2/pi]);