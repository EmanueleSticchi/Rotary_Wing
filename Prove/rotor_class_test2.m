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
rotore1 = rotore1.BEMT_articulated('T',W,V_inf,Chi,f);
W1 = W*1.3;
rotore1 = rotore1.BEMT_articulated('T',W1,V_inf,Chi,f);
W2 = W*1.5;
rotore1 = rotore1.BEMT_articulated('T',W2,V_inf,Chi,f);


%% Graphics
figure;
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{1,1}.Pc_Vec,'-k');
hold on
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{2,1}.Pc_Vec,':k');
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{3,1}.Pc_Vec,'.-k');

figure;
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{1,1}.beta0_Vec,'-k');
hold on
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{1,1}.beta1c_Vec,':k');
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{1,1}.beta1s_Vec,'.-k');
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{2,1}.beta0_Vec,'-r');
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{2,1}.beta1c_Vec,':r');
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{2,1}.beta1s_Vec,'.-r');
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{3,1}.beta0_Vec,'-g');
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{3,1}.beta1c_Vec,':g');
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{3,1}.beta1s_Vec,'.-g');

