%%  Sentiero di stallo per il rotore principale
clc; clear; close all
global aero
%% Data -------------------------------------------------------------------
rotore   = Rotor();
% aerodynamics
load('polari HH_02\polari\Aero_HH02_Re1250.mat')
rotore.Cl_alpha = aero.Cl_a_mode*180/pi;
rotore.Cl       = @(alpha) CL_(alpha);
rotore.Cd       = @(alpha) CD_(alpha);
alpha_max_2D = convang(15,'deg','rad');

% working conditions
rotore   = rotore.rot_vel('omega',726/24);  % [rad/s]
rotore.h = 0;                        % Assume density air = 1.23 Kg/m^3
rotore   = rotore.ambient();         % compute ambient conditions
V_inf    = linspace(0.1,convvel(500,'km/h','m/s'));
Chi      = convang(5,'deg','rad');
f        = 2;
W        = convforce(17650,'lbf','N');

% geometry
rotore   = rotore.r(linspace(0.1,1,50));      % dominio radiale [\]
rotore.R = convlength(24,'ft','m');           % Raggio rotore   [m]
rotore.N = 4;                                 % numero di pale  [\]
rotore.c = linspace(0.53,0.53,rotore.n_r);    % corde           [m]
I_MR     = 3800;
I_MR     = convmass(I_MR,'slug','kg');
I_MR     = convlength(convlength(I_MR,'ft','m'),'ft','m');
rotore   = rotore.mass_prop('I',I_MR);        % Mom. di inerzia [Kg*m^2]
rotore.theta_t = convang(-9,'deg','rad');     % theta twist     [deg]


[s,r,c] = rotore.sentiero_stallo(alpha_max_2D,1:0.01:1.12,'T',W,Chi,f);






%%  Graphics --------------------------------------------------------------
M_e = s.Mach_e;
figure
% Create polar data
[r,psi] = meshgrid(rotore.r_bar,s.options.Psi);
% Convert to Cartesian
x = r.*cos(psi);
y = r.*sin(psi);
% define polar axes
h = polar(x,y);
hold on;
polar(s.options.Psi,rotore.r_bar(1)*ones(length(s.options.Psi),1),'k')
polar(s.options.Psi,rotore.r_bar(end)*ones(length(s.options.Psi),1),'k')
% contourf(x,y,alpha_e');
pc= pcolor(x,y,M_e');
contour(x,y,M_e','k','ShowText','on');
shading interp
% colormap 'hsv'
cbar=colorbar(gca);
cbar.Label.String = 'M_e';
cbar.Label.FontSize= 16;
% cbar.Limits = [-10 10];

% Hide the POLAR function data and leave annotations
set(h,'Visible','off')
% Turn off axes and set square aspect ratio
axis off
axis image
view([90 90])
s_mu = sprintf('%0.2f',s.mu/cos(s.alpha_TPP_Vec));
s_a  = sprintf('%0.2f',max(M_e,[],'all'));
title(['\mu = ',s_mu,'   M_{e_{max}} = ',s_a])

