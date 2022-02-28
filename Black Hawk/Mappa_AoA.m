%%  Mappa per gli angoli per il rotore principale
clc; clear; close all
global aero
%% Data -------------------------------------------------------------------
rotore   = Rotor();
% aerodynamics
load('polari HH_02\polari\Aero_HH02_Re1250.mat')
rotore.Cl_alpha = aero.Cl_a_mode;
rotore.Cl       = @(alpha) CL_(alpha);
rotore.Cd       = @(alpha) CD_(alpha);

% working conditions
rotore   = rotore.rot_vel('omega',convvel(726,'ft/s','m/s')/14.63);
rotore.h = 0;                        % Assume density air = 1.23 Kg/m^3
rotore   = rotore.ambient();         % compute ambient conditions
V_inf    = linspace(0.1,convvel(200,'km/h','m/s'));
Chi      = convang(0.001,'deg','rad');
f        = 2;
W        = convforce(3000,'lbf','N');

% geometry
rotore   = rotore.r(linspace(0.1,1,20));      % dominio radiale [\]
rotore.R = 14.63;                             % Raggio rotore   [m]
rotore.N = 4;                                 % numero di pale  [\]
rotore.c = linspace(0.53,0.53,rotore.n_r);    % corde           [m]
I_MR     = 3800;
I_MR     = convmass(I_MR,'slug','kg');
I_MR     = convlength(convlength(I_MR,'ft','m'),'ft','m');
rotore   = rotore.mass_prop('I',I_MR);        % Mom. di inerzia [Kg*m^2]
rotore.theta_t = convang(-9,'deg','rad');     % theta twist     [deg]


%%  Graphics --------------------------------------------------------------

alphamap(rotore,'Solve',{'T',W,V_inf([1:10:100]),Chi,f,BEMTset_rotor()});

