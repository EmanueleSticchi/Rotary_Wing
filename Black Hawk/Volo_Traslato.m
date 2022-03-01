%%  Analisi Rotore Principale in volo tralato
clc; clear; close all
global aero
%% Data -------------------------------------------------------------------
rotore   = Rotor();
% aerodynamics
load('polari HH_02\polari\Aero_HH02_Re1250.mat')
rotore.Cl_alpha = aero.Cl_a_mode*180/pi;
rotore.Cl       = @(alpha) CL_(alpha);
rotore.Cd       = @(alpha) CD_(alpha);

% working conditions
rotore   = rotore.rot_vel('omega',726/24); % [rad/s]
rotore.h = 0;                        % Assume density air = 1.23 Kg/m^3
rotore   = rotore.ambient();         % compute ambient conditions
V_inf    = linspace(0.1,convvel(200,'km/h','m/s'));
Chi      = convang(0.001,'deg','rad');
f        = 2;
W        = convforce(3000,'lbf','N');

% geometry
rotore   = rotore.r(linspace(0.1,1,20));      % dominio radiale [\]
rotore.R = convlength(24,'ft','m');           % Raggio rotore   [m]
rotore.N = 4;                                 % numero di pale  [\]
rotore.c = linspace(0.53,0.53,rotore.n_r);    % corde           [m]
I_MR     = 3800;
I_MR     = convmass(I_MR,'slug','kg');
I_MR     = convlength(convlength(I_MR,'ft','m'),'ft','m');
rotore   = rotore.mass_prop('I',I_MR);        % Mom. di inerzia [Kg*m^2]
% rotore   = rotore.mass_prop('G',7);
rotore.theta_t = convang(-9,'deg','rad');     % theta twist     [deg]


%%  Analisys --------------------------------------------------------------
rotore = rotore.BEMT_articulated('T',W,V_inf,Chi,f,BEMTset_rotor());

%% Graphics ---------------------------------------------------------------
s = rotore.Analisi_articulated{1,1};
figure;
plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.Pc_Vec,'-k');
hold on

figure;
plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.Tc,'-k');
hold on

figure;
plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.theta0*180/pi,'-k');
hold on

figure;
plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.alpha_TPP_Vec*180/pi,'-k');
hold on

figure;
plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.lam_Vec,'-k');
hold on
plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.lam_i,':k');
plot(s.mu,s.lam_c,'.-k');

figure;
plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.beta0_Vec*180/pi,'-k');
hold on
plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.beta1c_Vec*180/pi,':k');
plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.beta1s_Vec*180/pi,'.-k');




