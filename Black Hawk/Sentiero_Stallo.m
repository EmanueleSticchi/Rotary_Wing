%%  Sentiero di stallo per il rotore principale
clc; clear; close all
global aero
pngflag = 1;
folder = 'Immagini\SentieroDiStallo\';
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


% [s,r,c] = rotore.sentiero_stallo(alpha_max_2D,1:0.01:1.12,'T',W,Chi,f);

VVec  = linspace(1,1.12,6);
counter = 0;

for i = 1:length(VVec)
    [s,r,c] = rotore.sentiero_stallo(alpha_max_2D,VVec(1:i),'T',W,Chi,f);
end

if pngflag == 1
    for i = 1:2*length(VVec)
        figure(i)
        if mod(i,2) == 0 || i == 1
            FileName = [folder,sprintf('Sds%d.eps', counter)];
            ax = gca;
            exportgraphics(ax,FileName)
            counter = counter + 1;
        end
    end
    else
        disp('Non stai generando nessun file .png!');
end

rotore.alphamap('Plot',{s;s.mu},'no')
ax = gca;
exportgraphics(ax,[folder,'Sds0.eps']);
