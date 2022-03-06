%%  Analisi di un'elica con la BEMT
clc; clear; close all
global aero
%%  DATI ----------------------------------------------------------
toll=1e-6;
el=Elica();

% Geometrici
% interpolazione dei dati effettivi
r    = convlength([10,12,18,24,30,36,42,48,54,60],'in','m');
rb   = r/r(end);
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

% figure
% plot(el.r_bar,el.c);
% figure
% plot(el.r_bar,el.theta*180/pi);

% Di funzionamento
theta_75 = convang(35,'deg','rad');
[~,idx] = min(abs(el.r_bar - 0.75)); 
el.theta = el.theta - el.theta(idx) + theta_75;  % setto il calettamento nominale al 75 %
el=el.rot_vel('RPM',1200);   % valore medio di quelli nel paper
el=el.altitude(0);
J=[0.1:0.05:1.5];
% Aerodinamici
load Aero_NACA16212.mat
el.Cl=@(alpha,r_bar,M,Re) CL_(alpha);
el.Cd=@(alpha,r_bar,M,Re) CD_(alpha);
% el.Cl = @(alpha,r_bar,M,Re) 2*pi*alpha;
% el.Cd = @(alpha,r_bar,M,Re) 0.02*alpha./alpha;

%% Analisi BEMT --------------------------------------------------------
alpha0=-10*pi/180;
alpha1=20*pi/180;
options=BEMTset();
% el=el.BEMT(J,alpha0,alpha1,options);
options.P_correction='on';
el=el.BEMT(J,alpha0,alpha1,options);
% el.LAMBDA(el.r_bar > 0.9 ) = linspace(0,2.5*pi/180,sum(el.r_bar >0.9)); 
% el.LAMBDA = linspace(0,20*pi/180,el.n_r);
% el=el.BEMT(J,alpha0,alpha1,options);
%% Post - Processing
an=1;
for an=1:1
figure(1)
plotta(J,el.Analisi{an, 1}.CT,{'J = $ \frac{V_{\infty}}{nD}$';'$C_T$ = $\frac{T}{\rho n^2 D^4}$'})
yline(0);

figure(2)
plotta(J,el.Analisi{an, 1}.CP,{'J = $ \frac{V_{\infty}}{nD}$';'$C_P$ = $\frac{P}{\rho n^3 D^5}$'})
yline(0);

figure(3)
log=el.Analisi{an, 1}.CT>=0;
plotta(J(log),el.Analisi{1, 1}.eta(log,1),{'J = $ \frac{V_{\infty}}{nD}$';'$\eta$ = $\frac{TV_{\infty}}{P}$'})
yline(0);
%
idx=10;
figure(4)
plotta(el.r_bar,el.Analisi{an,1}.dCt_dr_bar(idx,:),{'$ \bar{r}$';'$\frac{dC_T}{d\bar{r}}$'})
yline(0);
figure(5)
plotta(el.r_bar,el.Analisi{an,1}.dCq_dr_bar(idx,:),{'$ \bar{r}$';'$\frac{dC_Q}{d\bar{r}}$'})
yline(0);
figure(6)
plotta(el.r_bar,el.Analisi{an,1}.dCp_dr_bar(idx,:),{'$ \bar{r}$';'$\frac{dC_P}{d\bar{r}}$'})
yline(0);
figure(7)
plotta(el.r_bar,el.Analisi{an,1}.eta_e(idx,:),{'$ \bar{r}$';'$\eta_e$'})
yline(0);

figure(8)
plotta(el.r_bar,el.Analisi{an,1}.Mach(idx,:),{'$ \bar{r}$';'$M$'})
figure(9)
plotta(el.r_bar,el.Analisi{an,1}.Re(idx,:),{'$ \bar{r}$';'$Re$'})
end

%% MODEL 3D
data=importdata('NACA 16-212.dat');
x=data.data(:,1);
z=data.data(:,2);
figure
el.Model3D(x,z)


%% Function ---------------------------------------------------------
function CL=CL_(alpha)
global aero

    alpha=alpha*180/pi;
    if alpha <  aero.alpha(1)
        CL=aero.Cl(1);
    elseif alpha > aero.alpha(end)
        CL=aero.Cl(end);
    else 
        CL=interp1(aero.alpha,aero.Cl,alpha);
    end
end

function CD=CD_(alpha)
global aero
    alpha=alpha*180/pi;
    if alpha <  aero.alpha(1)
        CD=aero.Cd(1);
    elseif alpha > aero.alpha(end)
        CD=aero.Cd(end);
    else 
        CD=interp1(aero.alpha,aero.Cd,alpha);
    end
end