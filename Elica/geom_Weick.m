%%  Geometry Weick
clc; clear; close all
folder = 'immagini/geom/';
%%  DATI ----------------------------------------------------------
toll=1e-6;
el=Elica();

% Geometrici
% interpolazione dei dati effettivi
r    = convlength([10,12,18,24,30,36,42,48,54,60],'in','m');
rb   = r/convlength(63,'in','m');
c    = convlength([6.97,...     % Corda delle sezioni, [m]
    7.36,8.08,8.41,8.43,8.17,7.53,6.53,5.21,3.74],'in','m');
theta = pi/180*....             % Angolo di calettamento, [rad]
    [42.4,39.1,35.1,30.2,26.7,23.9,21.75,20.1,18.8,17.3];
c_   = polyfit(rb,c,5);         % definizione di una polinomiale di 5*grado
t_   = polyfit(rb,theta,5);     % definizione di una polinomiale di 5*grado

el   = el.r_(rb(1):0.01:1);     % Definizione dominio radiale
el.N = 3;                       % Numero di pale    
el.R = convlength(63,'in','m'); % Raggio dell'elica, [m]
el.c = polyval(c_,el.r_bar);    % Corda delle sezioni, [m]
el   = el.derived_properties;   % compute some derived prop. like as sigma

el.theta = polyval(t_,el.r_bar); % Angolo di calettamento, [rad]
el.LAMBDA = zeros(el.n_r,1);    % Angolo di freccia, [rad]

% Di funzionamento
theta_75 = convang(0,'deg','rad');
[~,idx] = min(abs(el.r_bar - 0.75)); 
el.theta = el.theta - el.theta(idx) + theta_75;  % setto il calettamento nominale al 75 %
el=el.rot_vel('RPM',1200);   % valore medio di quelli nel paper
el=el.altitude(0);
J=[0.05:0.05:1.3];
% Aerodinamici
el.Cl = @(alpha,r_bar,M,Re) 2*pi*alpha;
el.Cd = @(alpha,r_bar,M,Re) 0.02*alpha./alpha;
%% MODEL 3D
data=importdata('NACA 16-212.dat');
x=data.data(:,1);
z=data.data(:,2);
el.Model3D(x,z)
el.Model3D(x,z)
figure(1)
ax = gca;
view(ax,[-90 90]); 
f.Units = 'normalized';
f.Position = [0.13,0.11,0.775,0.815];
figure(2)
ax = gca;
view(ax,[25 30]); 
axis off

figure
plotta(el.r_bar,el.c,{'$ \bar{r}$';'c [m]'});
figure
plotta(el.r_bar,el.theta*180/pi,{'$ \bar{r}$';'$\theta$ [deg]'});

%% Save
count = 0;
for i =1:6
    count = count + 1;
    figure(i)
    FileName = sprintf(['geom','%d.eps'], count);
    ax = gca;
    exportgraphics(ax,[folder,FileName])
end
