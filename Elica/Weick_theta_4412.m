%%  Analisi di un'elica Weick con la BEMT al varia del passo di radice
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

% Aerodinamici
load Aero_NACA4412.mat
el.Cl=@(alpha,r_bar,M,Re) CL_(alpha);
el.Cd=@(alpha,r_bar,M,Re) CD_(alpha);

% Di funzionamento
vtheta = convang(linspace(15,50,5),'deg','rad');
el=el.rot_vel('RPM',1200);   % valore medio di quelli nel paper
el=el.altitude(0);
options=BEMTset();      options.P_correction='on'; 
% alpha0 = [-2 -10 -10 -30 -50]*pi/180;
alpha0 = -2*pi/180;
alpha1 = [10 10 20 20 20]*pi/180;
J_end  = [1.3 1.4 1.6 2 2.5];
for tdx = 1:length(vtheta)
    theta_75 = vtheta(tdx);
    [~,idx] = min(abs(el.r_bar - 0.75));
    el.theta = el.theta - el.theta(idx) + theta_75;  % setto il calettamento nominale al 75 %
    J=[0.1:0.05:J_end(tdx)];
    %% Analisi BEMT --------------------------------------------------------       
    el=el.BEMT(J,alpha0,alpha1(tdx),options);
    
end

%% Post - Processing    
formatspec = {'-k';'--k';':k';'.-k';'^-k';'s-k'};
for tdx = 1:length(vtheta)
    J=[0.1:0.05:J_end(tdx)];
    s = el.Analisi{tdx, 1};
    name = ['\theta_{75} = ',num2str(vtheta(tdx)*180/pi)];
    figure(1)
    plotta(J,s.CT,...
        {'J = $ \frac{V_{\infty}}{nD}$';'$C_T$ = $\frac{T}{\rho n^2 D^4}$'},...
        formatspec{tdx,1},name)
    
    figure(2)
    plotta(J,s.CP,...
        {'J = $ \frac{V_{\infty}}{nD}$';'$C_P$ = $\frac{P}{\rho n^3 D^5}$'},...
        formatspec{tdx,1},name)
    
    figure(3)
    log=s.CT>=0;
    plotta(J(log),s.eta(log,1),...
        {'J = $ \frac{V_{\infty}}{nD}$';'$\eta$ = $\frac{TV_{\infty}}{P}$'},...
        formatspec{tdx,1},name)    
end

for i =1:3
    figure(i)
    lg = legend();
    lg.AutoUpdate='off';
    yline(0);
    xline(0);
end


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