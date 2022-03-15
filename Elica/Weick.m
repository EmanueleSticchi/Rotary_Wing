%%  Analisi di un'elica Weick con la BEMT
clc; clear; close all
load CC_experimential_Weick.mat
flag = 1;
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
theta_75 = convang(15.5,'deg','rad');
[~,idx] = min(abs(el.r_bar - 0.75)); 
el.theta = el.theta - el.theta(idx) + theta_75;  % setto il calettamento nominale al 75 %
el=el.rot_vel('RPM',1200);   % valore medio di quelli nel paper
el=el.altitude(0);
J=[0.02:0.05:1];
% Aerodinamici
el.Cl = @(alpha,r_bar,M,Re) 2*pi*alpha;
el.Cd = @(alpha,r_bar,M,Re) 0.01;

%% Analisi BEMT --------------------------------------------------------

alpha0=-2*pi/180;
options=BEMTset();      options.P_correction='on';
options.Hub_correction = 'on'; options.Cd_hub =0.5;
el=el.BEMT(J,alpha0,options);

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
    lg.Color = 'none';
    
    if sum([4,5,6,12] == i)
        lg.Location = 'northwest';
    elseif i==11
        lg.Location = 'southwest';
    else
        lg.Location = 'southeast';
    end
    yline(0);
    xline(0);
end


%% Salvataggio

if flag ==1
    save(1:3 ,'car','immagini/Weick_2pi/caratteristiche/');
    save(4:15,'distr','immagini/Weick_2pi/distribuzioni/');
end

%% Plot comparison with experimental data
n_num = 'Numerico';     n_exp = 'Sperimentale';
figure
plotta(J,s.CT,...
    {'J = $ \frac{V_{\infty}}{nD}$';'$C_T$ = $\frac{T}{\rho n^2 D^4}$'},...
    '-k',n_num)
plotta(CC_exp.J,CC_exp.Ct,...
    {'J = $ \frac{V_{\infty}}{nD}$';'$C_T$ = $\frac{T}{\rho n^2 D^4}$'},...
    '--k',n_exp)


figure
plotta(J,s.CP,...
    {'J = $ \frac{V_{\infty}}{nD}$';'$C_P$ = $\frac{P}{\rho n^3 D^5}$'},...
    '-k',n_num)
plotta(CC_exp.J,CC_exp.Cp,...
    {'J = $ \frac{V_{\infty}}{nD}$';'$C_P$ = $\frac{P}{\rho n^3 D^5}$'},...
    '--k',n_exp)


figure
log=s.CT>=0;
plotta(J(log),s.eta(log,1),...
    {'J = $ \frac{V_{\infty}}{nD}$';'$\eta$ = $\frac{TV_{\infty}}{P}$'},...
    '-k',n_num)
plotta(CC_exp.J,CC_exp.eta,...
    {'J = $ \frac{V_{\infty}}{nD}$';'$\eta$ = $\frac{TV_{\infty}}{P}$'},...
    '--k',n_exp)

for i =16:18
    figure(i)
    lg = legend();
    lg.AutoUpdate='off';
    lg.Color = 'none';
    lg.Location = 'best';
    yline(0);
    xline(0);
end

if flag ==1
    save(16:18 ,'car_comp','immagini/Weick_2pi/caratteristiche/');
end



% % %% MODEL 3D
% data=importdata('NACA 16-212.dat');
% x=data.data(:,1);
% z=data.data(:,2);
% el.Model3D(x,z)
% 
% figure
% plotta(el.r_bar,el.c,{'$ \bar{r}$';'c [m]'});
% figure
% plotta(el.r_bar,el.theta*180/pi,{'$ \bar{r}$';'$\theta$ [deg]'});


%% Function 

function save(idxF,prename,folder)
count = 0;
for i =1:length(idxF)
    count = count + 1;
    figure(idxF(i))
    %     legend
    FileName = sprintf([prename,'%d.eps'], count);
    ax = gca;
    exportgraphics(ax,[folder,FileName])
end
end