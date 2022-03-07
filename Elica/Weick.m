%%  Analisi di un'elica Weick con la BEMT
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

% Di funzionamento
theta_75 = convang(15,'deg','rad');
[~,idx] = min(abs(el.r_bar - 0.75)); 
el.theta = el.theta - el.theta(idx) + theta_75;  % setto il calettamento nominale al 75 %
el=el.rot_vel('RPM',1200);   % valore medio di quelli nel paper
el=el.altitude(0);
J=[0.05:0.05:1.3];
% Aerodinamici
el.Cl = @(alpha,r_bar,M,Re) 2*pi*alpha;
el.Cd = @(alpha,r_bar,M,Re) 0.01*alpha./alpha;

%% Analisi BEMT --------------------------------------------------------

alpha0=2*pi/180;      alpha1=10*pi/180;
options=BEMTset();      options.P_correction='on';
el=el.BEMT(J,alpha0,alpha1,options);

%% Post - Processing
s = el.Analisi{1, 1};
figure(1)
plotta(J,s.CT,{'J = $ \frac{V_{\infty}}{nD}$';'$C_T$ = $\frac{T}{\rho n^2 D^4}$'})
yline(0);

figure(2)
plotta(J,s.CP,{'J = $ \frac{V_{\infty}}{nD}$';'$C_P$ = $\frac{P}{\rho n^3 D^5}$'})
yline(0);

figure(3)
log=s.CT>=0;
plotta(J(log),s.eta(log,1),{'J = $ \frac{V_{\infty}}{nD}$';'$\eta$ = $\frac{TV_{\infty}}{P}$'})
yline(0);

%
formatspec = {'-k';'--k';':k';'.-k';'^-k';'s-k'}; k = 0;
for jdx = 1:length(J)
    
    if mod(jdx,5) == 0
        k=k+1;
        name = ['J = ',num2str(J(jdx))];
        figure(4)
        plotta(el.r_bar,s.dCt_dr_bar(jdx,:),...
            {'$ \bar{r}$';'$\frac{dC_T}{d\bar{r}}$'},...
            formatspec{k,1},name)

        figure(5)
        plotta(el.r_bar,s.dCq_dr_bar(jdx,:),...
            {'$ \bar{r}$';'$\frac{dC_Q}{d\bar{r}}$'},...
            formatspec{k,1},name)

        figure(6)
        plotta(el.r_bar,s.dCp_dr_bar(jdx,:),...
            {'$ \bar{r}$';'$\frac{dC_P}{d\bar{r}}$'},...
            formatspec{k,1},name)

        figure(7)
        plotta(el.r_bar,s.eta_e(jdx,:),{'$ \bar{r}$';'$\eta_e$'},...
            formatspec{k,1},name)


        figure(8)
        plotta(el.r_bar,s.Mach(jdx,:),{'$ \bar{r}$';'$M$'},...
            formatspec{k,1},name)
        figure(9)
        plotta(el.r_bar,s.Re(jdx,:),{'$ \bar{r}$';'$Re$'},...
            formatspec{k,1},name)
        figure(10)
        plotta(el.r_bar,s.alpha(jdx,:)*180/pi,{'$ \bar{r}$';'$\alpha$'},...
            formatspec{k,1},name)
        figure(11)
        plotta(el.r_bar,s.phi(jdx,:)*180/pi,{'$ \bar{r}$';'$\varphi$'},...
            formatspec{k,1},name)
        figure(12)
        plotta(el.r_bar,s.a(jdx,:),{'$ \bar{r}$';'$a$'},...
            formatspec{k,1},name)
        figure(13)
        plotta(el.r_bar,s.ap(jdx,:),{'$ \bar{r}$';'$a$'''},...
            formatspec{k,1},name)
        figure(14)
        plotta(el.r_bar,s.lambda1(jdx,:),{'$ \bar{r}$';'$\lambda_1$'},...
            formatspec{k,1},name)
        figure(15)
        plotta(el.r_bar,s.lambda2(jdx,:),{'$ \bar{r}$';'$\lambda_2$'},...
            formatspec{k,1},name)
    end
end
for i =4:15
    figure(i)
    lg = legend();
    lg.AutoUpdate='off';
    yline(0);
    xline(0);
end


%% MODEL 3D
data=importdata('NACA 16-212.dat');
x=data.data(:,1);
z=data.data(:,2);
figure
el.Model3D(x,z)

figure
plotta(el.r_bar,el.c,{'$ \bar{r}$';'c [m]'});
figure
plotta(el.r_bar,el.theta*180/pi,{'$ \bar{r}$';'$\theta$ [deg]'});