%% Test BEMT SALITA 
% Confronto con Eurocopter AS365N fatto in classe su Excel
% Esito positivo.
clear all; close all; clc;
%% dati
rotore1 = Rotor();
% working conditions and other inputs
dim_vel = 50;
V_inf   = linspace(0,convvel(10,'km/h','m/s'),dim_vel); % climb speed
rotore1.h     = 0;
rotore1 = rotore1.ambient();

% properties
rotore1 = rotore1.r(linspace(0.1,1,20));
rotore1.R     = 11.94/2;
rotore1.N     = 4;
rotore1.c     = linspace(0.39,0.39,rotore1.n_r);
theta0 = 16;    Dtheta=10;
rotore1.theta = pi/180*(theta0 - Dtheta*rotore1.r_bar);
rotore1 = rotore1.derived_properties();
rotore1 = rotore1.rot_vel('omega',238/rotore1.R);
rotore1.Cl = @(alpha) 5.7*alpha;

rotore1 = rotore1.BEMT_salita(V_inf);


%% Graphics
n=1;
s=rotore1.Analisi_salita{n,1};

figure;
plot(s.mu,s.Tc)
xlabel('\mu')
ylabel('T_c')

figure;
plot(s.mu,s.Qc)
xlabel('\mu')
ylabel('Q_c')

idx_V = 1;
figure;
plot(rotore1.r_bar,s.dTc(idx_V,:))
xlabel('\bar{r}')
ylabel('dT_c/dr')

figure;
plot(rotore1.r_bar,s.dQc(idx_V,:))
xlabel('\bar{r}')
ylabel('dQ_c/dr')


