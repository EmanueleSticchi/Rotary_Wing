clear all; close all; clc;

rotore1 = Rotor();
% working conditions and other inputs
dim_vel = 50;
V_inf   = linspace(0.1,convvel(293,'km/h','m/s'),dim_vel);
Chi     = convang(5,'deg','rad');
f       = 3;
W       = 75278;
rotore1.theta_t = convang(-8,'deg','rad');
rotore1.h     = 0;
% properties
rotore1 = rotore1.r(linspace(0.1,1,20));
rotore1.R     = 7.6;
rotore1.N     = 3;
rotore1.c     = linspace(0.4,0.4,rotore1.n_r);
% rotore1.theta = pi/180*linspace(13.3,9,rotore1.n_r);
% obj.I = convmass(obj.Ixx,'slug','kg');
% obj.I = convlength(convlength(obj.Ixx,'ft','m'),'ft','m');
% function recall
rotore1 = rotore1.ambient();

rotore1 = rotore1.mass_prop('G',8);
rotore1 = rotore1.rot_vel('omega',1,1);
theta0  = convang(15,'deg','rad');  
rotore1 = rotore1.BEMT_articulated('Theta',theta0,V_inf,Chi,f);
theta0  = convang(17,'deg','rad');
rotore1 = rotore1.BEMT_articulated('Theta',theta0,V_inf,Chi,f);
theta0  = convang(20,'deg','rad');
rotore1 = rotore1.BEMT_articulated('Theta',theta0,V_inf,Chi,f);


%% Graphics
figure;
plot(rotore1.Analisi_articulated{1,1}.mu,rotore1.Analisi_articulated{1,1}.Pc_Vec,'-k');
hold on
plot(rotore1.Analisi_articulated{2,1}.mu,rotore1.Analisi_articulated{2,1}.Pc_Vec,':k');
plot(rotore1.Analisi_articulated{3,1}.mu,rotore1.Analisi_articulated{3,1}.Pc_Vec,'.-k');

figure;
plot(rotore1.Analisi_articulated{1,1}.mu,rotore1.Analisi_articulated{1,1}.beta0_Vec,'-k');
hold on
plot(rotore1.Analisi_articulated{1,1}.mu,rotore1.Analisi_articulated{1,1}.beta1c_Vec,':k');
plot(rotore1.Analisi_articulated{1,1}.mu,rotore1.Analisi_articulated{1,1}.beta1s_Vec,'.-k');
plot(rotore1.Analisi_articulated{2,1}.mu,rotore1.Analisi_articulated{2,1}.beta0_Vec,'-r');
plot(rotore1.Analisi_articulated{2,1}.mu,rotore1.Analisi_articulated{2,1}.beta1c_Vec,':r');
plot(rotore1.Analisi_articulated{2,1}.mu,rotore1.Analisi_articulated{2,1}.beta1s_Vec,'.-r');
plot(rotore1.Analisi_articulated{3,1}.mu,rotore1.Analisi_articulated{3,1}.beta0_Vec,'-g');
plot(rotore1.Analisi_articulated{3,1}.mu,rotore1.Analisi_articulated{3,1}.beta1c_Vec,':g');
plot(rotore1.Analisi_articulated{3,1}.mu,rotore1.Analisi_articulated{3,1}.beta1s_Vec,'.-g');

