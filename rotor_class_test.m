clear all; close all; clc;

rotore1 = Rotor();
% working conditions
dim_vel = 50;
V_inf   = linspace(0,200,dim_vel);

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
options = BEMTset_rotor();
rotore1 = rotore1.BEMT_rotore(V_inf,options);

figure;
for i = 1:dim_vel
    plot(rotore1.r_bar,rotore1.Analisi{1,1}.dTc(i,:),'-k')
    hold on
end
figure;
plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi{1,1}.Tc,'-k');

