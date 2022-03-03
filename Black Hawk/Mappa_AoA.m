%%  Mappa per gli angoli per il rotore principale
clc; clear; close all
global aero
pngflag = 1;       %flag per il salvatagio delle immagini 
folderA = 'Immagini\Mappa_AoA\';
folderM = 'Immagini\Mappa_Mach\';
%% Data -------------------------------------------------------------------
rotore   = Rotor();
% aerodynamics
load('polari HH_02\polari\Aero_HH02_Re1250.mat')
rotore.Cl_alpha = aero.Cl_a_mode*180/pi;
rotore.Cl       = @(alpha) CL_(alpha);
rotore.Cd       = @(alpha) CD_(alpha);

% working conditions
rotore   = rotore.rot_vel('omega',726/24);  % [rad/s]
rotore.h = 0;                        % Assume density air = 1.23 Kg/m^3
rotore   = rotore.ambient();         % compute ambient conditions
V_inf    = linspace(0.1,convvel(350,'km/h','m/s'),6);
Chi      = convang(5,'deg','rad');
f        = 2;
W        = convforce(17650,'lbf','N');

% geometry
rotore   = rotore.r(linspace(0.1,1,20));      % dominio radiale [\]
rotore.R = convlength(24,'ft','m');           % Raggio rotore   [m]
rotore.N = 4;                                 % numero di pale  [\]
rotore.c = linspace(0.53,0.53,rotore.n_r);    % corde           [m]
I_MR     = 3800;
I_MR     = convmass(I_MR,'slug','kg');
I_MR     = convlength(convlength(I_MR,'ft','m'),'ft','m');
rotore   = rotore.mass_prop('I',I_MR);        % Mom. di inerzia [Kg*m^2]
rotore.theta_t = convang(-9,'deg','rad');     % theta twist     [deg]


%%  Graphics --------------------------------------------------------------

alphamap(rotore,'Solve',{'T',W,V_inf,Chi,f,BEMTset_rotor()});
if pngflag == 1
    counter = 0;
    for i = 1:length(V_inf)
        counter = counter + 1;
        figure(i)
        FileName = sprintf('AoA%d.eps', counter);
        ax = gca;
        exportgraphics(ax,[folderA,FileName])
    end
else
    disp('Non stai generando nessun file .png!');
end





MachMap(rotore,'Solve',{'T',W,V_inf,Chi,f,BEMTset_rotor()});
% il numero di mach supera il valore sonico per velocità minore della V_ne.
% In teoria bisognerebbe tenere in conto dell'effetto della freccia
if pngflag == 1
    counter = 0;
    for i = (length(V_inf)+1):(2*length(V_inf))
        counter = counter + 1;
        figure(i)
        FileName = sprintf('Mach%d.eps', counter);
        ax = gca;
        exportgraphics(ax,[folderM,FileName])
    end
else
    disp('Non stai generando nessun file .png!');
end

