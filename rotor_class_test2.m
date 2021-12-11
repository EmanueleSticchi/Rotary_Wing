clear all; close all; clc;

rotore1 = Rotor();
% working conditions and other inputs
dim_vel = 50;
V_inf   = linspace(0,convvel(293,'km/h','m/s'),dim_vel);
Chi     = convang(5,'deg','rad');
f       = 3;
W       = 75278;
theta_t = convang(-8,'deg','rad');
rotore1.h     = 0;
% properties
rotore1 = rotore1.r(linspace(0.1,1,20));
rotore1.R     = 7.6;
rotore1.N     = 3;
rotore1.c     = linspace(0.4,0.4,rotore1.n_r);
rotore1.theta = pi/180*linspace(13.3,9,rotore1.n_r);
% function recall
rotore1 = rotore1.ambient();

rotore1 = rotore1.derived_properties();
rotore1 = rotore1.rot_vel('omega',1,1);
rotore1 = rotore1.articulated_rotor(V_inf,Chi,f,W,theta_t);
%% Graphics
figure;
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{1,1}.Pc_Vec.*10^3,'-k');
