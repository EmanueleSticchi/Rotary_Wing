%%  Mappa per gli angoli per il rotore principale
clc; clear; close all
global aero
pngflag = 0;       %flag per il salvatagio delle immagini 
folderA = 'Immagini\Mappa_AoA\';
folderM = 'Immagini\Mappa_Mach\';
%% Data -------------------------------------------------------------------
rotore   = Rotor();
% aerodynamics
load('polari HH_02\polari\Aero_HH02_Re1250.mat')
rotore.Cl_alpha = aero.Cl_a_mode*180/pi;
% rotore.Cl       = @(alpha) CL_(alpha);
% rotore.Cd       = @(alpha) CD_(alpha);

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
%                    Alpha
% a colori
s = alphamap(rotore,'Solve',{'T',W,V_inf,Chi,f,BEMTset_rotor()});
idx = 1:length(V_inf);
salva('AoA',folderA,idx,pngflag)

% non colori
alphamap(rotore,'Plot',{s;s.mu},'no');
idx = idx + length(V_inf);
salva('AoA_nc',folderA,idx,pngflag)


%                                 Mach
% colori
MachMap(rotore,'Plot',{s;s.mu});
idx = idx + length(V_inf);
salva('Mach',folderM,idx,pngflag)

% non colori
MachMap(rotore,'Plot',{s;s.mu},'no');
idx = idx + length(V_inf);
salva('Mach_nc',folderM,idx,pngflag)

                            % Mach con freccia
                            
Freccia = zeros(rotore.n_r,1);
Freccia(rotore.r_bar > 0.97) = convang(45,'deg','rad');
for i =1:rotore.n_r
    s.Mach_e(i,:,:) = s.Mach_e(i,:,:)*cos(Freccia(i));
end
% non a colori
MachMap(rotore,'Plot',{s;s.mu},'no');
idx = idx + length(V_inf);
salva('Mach_nc_F',folderM,idx,pngflag)

% colori
MachMap(rotore,'Plot',{s;s.mu});
idx = idx + length(V_inf);
salva('Mach_F',folderM,idx,pngflag)

%% PLot rotore con freccia 
rotore.LAMBDA = Freccia;
data=importdata('polari HH_02\HH_02.dat');
x=data.data(:,1);
z=data.data(:,2);
rotore.Model3D(x,z,20*pi/180)

idx = idx(end) + [1,2];
salva('rotore','Immagini\rotore\',idx,pngflag)


%% Function
function salva(prename,folder,idx,pngflag)
if pngflag == 1
    counter = 0;
    for i = idx(1):idx(end)
        counter = counter + 1;
        figure(i)
        FileName = sprintf([prename,'%d.eps'], counter);
        ax = gca;
        exportgraphics(ax,[folder,FileName])
    end
else
    disp('Non stai generando nessun file .png!');
end

end