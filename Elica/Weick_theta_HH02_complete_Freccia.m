%%  Analisi di un'elica Weick con la BEMT al varia del passo di radice
clc; clear; close all
global m_alpha m_Re m_Cl m_Cd
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

el.theta  = polyval(t_,el.r_bar); % Angolo di calettamento, [rad]
el.LAMBDA = linspace(0,10*pi/180,el.n_r);    % Angolo di freccia, [rad]

% Aerodinamici
load AeroHH02_complete.mat
el.Cl=@(alpha,r_bar,M,Re) Cl_(alpha,M,Re);
el.Cd=@(alpha,r_bar,M,Re) Cd_(alpha,M,Re);

% Di funzionamento
vtheta = convang(linspace(15,50,5),'deg','rad');
el=el.rot_vel('RPM',1200);   % valore medio di quelli nel paper
el=el.altitude(0);
options=BEMTset();      options.P_correction='on'; 
% alpha0 = [-2 -10 -10 -30 -50]*pi/180;
alpha0 = -2*pi/180;
J_end  = [1.3 1.4 1.6 2 2.5];
for tdx = 1:length(vtheta)
    theta_75 = vtheta(tdx);
    [~,idx] = min(abs(el.r_bar - 0.75));
    el.theta = el.theta - el.theta(idx) + theta_75;  % setto il calettamento nominale al 75 %
    J=[0.1:0.05:J_end(tdx)];
    %% Analisi BEMT --------------------------------------------------------       
    el=el.BEMT(J,alpha0,options);
    
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
function Cl = Cl_(alpha,M,Re)
% angoli in radianti
    global m_alpha m_Re m_Cl
    alpha     = convang(alpha,'rad','deg');
    alpha_lim = [m_alpha(1,1),m_alpha(1,end)];
    Re_lim    = [m_Re(1,1)   ,m_Re(end,1)];
    
    Cl = interp2(m_alpha,m_Re,m_Cl,...
        min(max(alpha_lim(1),alpha),alpha_lim(2)),...
        min(max(Re_lim(1),Re),Re_lim(2)));
    if M <1
        Cl = Cl.*sqrt(1-M.^2).^-1;
    else
        Cl = NaN;
    end
end

function Cd = Cd_(alpha,M,Re)
% angoli in radianti
    global m_alpha m_Re m_Cd
    alpha     = convang(alpha,'rad','deg');
    alpha_lim = [m_alpha(1,1),m_alpha(1,end)];
    Re_lim    = [m_Re(1,1)   ,m_Re(end,1)];
    
    Cd = interp2(m_alpha,m_Re,m_Cd,...
        min(max(alpha_lim(1),alpha),alpha_lim(2)),...
        min(max(Re_lim(1),Re),Re_lim(2)));
    if M <1
        k  = 2.5;
        Cd = Cd/sqrt(1-M^k);
    else
        Cd = NaN;
    end
end