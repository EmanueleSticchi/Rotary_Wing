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
theta_75 = convang(20,'deg','rad');
[~,idx] = min(abs(el.r_bar - 0.75)); 
el.theta = el.theta - el.theta(idx) + theta_75;  % setto il calettamento nominale al 75 %
el=el.rot_vel('RPM',1200);   % valore medio di quelli nel paper
el=el.altitude(0);
J=[0.05:0.05:1];

%% Analisi BEMT --------------------------------------------------------
formatspec = {'-k';'--k';':k';'.-k';'^-k';'s-k'};
for i = 1:3
    name = ['Case ',num2str(i)];
    if i==1
        el.LAMBDA = zeros(el.n_r,1);
    elseif i==2
        el.LAMBDA = linspace(0,35*pi/180,el.n_r);
    elseif i== 3
        el.LAMBDA = ones(el.n_r,1)*40*pi/180;
%     elseif i == 4
%         el.LAMBDA = zeros(el.n_r,1);
%         el.LAMBDA(el.r_bar > 0.9) = 10*pi/180;
    end
    alpha0 = -5*pi/180;
    options=BEMTset();      options.P_correction='on';
    el=el.BEMT(J,alpha0,options);
    %% Post - Processing
    s = el.Analisi{el.n_analisi, 1};
    figure(1)
    hold on 
    plotta(J,s.CT,...
        {'J = $ \frac{V_{\infty}}{nD}$';'$C_T$ = $\frac{T}{\rho n^2 D^4}$'},...
        formatspec{i,1},name)

    
    figure(2)
    hold on
    plotta(J,s.CP,...
        {'J = $ \frac{V_{\infty}}{nD}$';'$C_P$ = $\frac{P}{\rho n^3 D^5}$'},...
        formatspec{i,1},name)
    figure(3)
    hold on
    log=s.CT>=0;
    plotta(J(log),s.eta(log,1),...
        {'J = $ \frac{V_{\infty}}{nD}$';'$\eta$ = $\frac{TV_{\infty}}{P}$'},...
        formatspec{i,1},name)
    
    figure
    data=importdata('NACA 16-212.dat');
    x=data.data(:,1);
    z=data.data(:,2);
    el.Model3D(x,z)
end
for i =1:4
    figure(i)
    lg = legend();
    lg.AutoUpdate='off';
    yline(0);
    xline(0);
end


%% Post - Processing
% s = el.Analisi{1, 1};
% figure(1)
% plotta(J,s.CT,{'J = $ \frac{V_{\infty}}{nD}$';'$C_T$ = $\frac{T}{\rho n^2 D^4}$'})
% yline(0);
% 
% figure(2)
% plotta(J,s.CP,{'J = $ \frac{V_{\infty}}{nD}$';'$C_P$ = $\frac{P}{\rho n^3 D^5}$'})
% yline(0);
% 
% figure(3)
% log=s.CT>=0;
% plotta(J(log),s.eta(log,1),{'J = $ \frac{V_{\infty}}{nD}$';'$\eta$ = $\frac{TV_{\infty}}{P}$'})
% yline(0);
% 
% 
% formatspec = {'-k';'--k';':k';'.-k';'^-k';'s-k'}; k = 0;
% for jdx = 1:length(J)
%     
%     if mod(jdx,5) == 0
%         k=k+1;
%         name = ['J = ',num2str(J(jdx))];
%         figure(4)
%         plotta(el.r_bar,s.dCt_dr_bar(jdx,:),...
%             {'$ \bar{r}$';'$\frac{dC_T}{d\bar{r}}$'},...
%             formatspec{k,1},name)
% 
%         figure(5)
%         plotta(el.r_bar,s.dCq_dr_bar(jdx,:),...
%             {'$ \bar{r}$';'$\frac{dC_Q}{d\bar{r}}$'},...
%             formatspec{k,1},name)
% 
%         figure(6)
%         plotta(el.r_bar,s.dCp_dr_bar(jdx,:),...
%             {'$ \bar{r}$';'$\frac{dC_P}{d\bar{r}}$'},...
%             formatspec{k,1},name)
% 
%         figure(7)
%         plotta(el.r_bar,s.eta_e(jdx,:),{'$ \bar{r}$';'$\eta_e$'},...
%             formatspec{k,1},name)
% 
% 
%         figure(8)
%         plotta(el.r_bar,s.Mach(jdx,:),{'$ \bar{r}$';'$M$'},...
%             formatspec{k,1},name)
%         figure(9)
%         plotta(el.r_bar,s.Re(jdx,:),{'$ \bar{r}$';'$Re$'},...
%             formatspec{k,1},name)
%         figure(10)
%         plotta(el.r_bar,s.alpha(jdx,:)*180/pi,{'$ \bar{r}$';'$\alpha$'},...
%             formatspec{k,1},name)
%         figure(11)
%         plotta(el.r_bar,s.phi(jdx,:)*180/pi,{'$ \bar{r}$';'$\varphi$'},...
%             formatspec{k,1},name)
%         figure(12)
%         plotta(el.r_bar,s.a(jdx,:),{'$ \bar{r}$';'$a$'},...
%             formatspec{k,1},name)
%         figure(13)
%         plotta(el.r_bar,s.ap(jdx,:),{'$ \bar{r}$';'$a$'''},...
%             formatspec{k,1},name)
%         figure(14)
%         plotta(el.r_bar,s.lambda1(jdx,:),{'$ \bar{r}$';'$\lambda_1$'},...
%             formatspec{k,1},name)
%         figure(15)
%         plotta(el.r_bar,s.lambda2(jdx,:),{'$ \bar{r}$';'$\lambda_2$'},...
%             formatspec{k,1},name)
%     end
% end
% for i =4:15
%     figure(i)
%     lg = legend();
%     lg.AutoUpdate='off';
%     yline(0);
%     xline(0);
% end


%% MODEL 3D
% data=importdata('NACA 16-212.dat');
% x=data.data(:,1);
% z=data.data(:,2);
% el.Model3D(x,z)
% 
% figure
% plotta(el.r_bar,el.c,{'$ \bar{r}$';'c [m]'});
% figure
% plotta(el.r_bar,el.theta*180/pi,{'$ \bar{r}$';'$\theta$ [deg]'});


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