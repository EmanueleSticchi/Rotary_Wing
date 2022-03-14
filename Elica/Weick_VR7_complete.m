%%  Analisi di un'elica Weick con la BEMT
% modello aerodinamico completo profilo BV VR7
clc; clear; close all
flag = 0;
global m_alpha m_Re m_Cl m_Cd
load CC_experimential_Weick.mat
%%  DATI ----------------------------------------------------------
% importo l'oggeto elica
load el_Weick

% Di funzionamento
theta_75 = convang(15.5,'deg','rad');
[~,idx] = min(abs(el.r_bar - 0.75)); 
el.theta = el.theta - el.theta(idx) + theta_75;  % setto il calettamento nominale al 75 %
el=el.rot_vel('RPM',1200);   % valore medio di quelli nel paper
el=el.altitude(0);
J=[0.05:0.05:0.9];
% Aerodinamici
load AeroVR7_complete.mat
el.Cl=@(alpha,r_bar,M,Re) Cl_(alpha,M,Re);
el.Cd=@(alpha,r_bar,M,Re) Cd_(alpha,M,Re);

%% Analisi BEMT --------------------------------------------------------

alpha0=-2*pi/180;
options=BEMTset();              options.P_correction = 'on';
options.Hub_correction = 'on';  options.Cd_hub       = 0.5;
options.toll = 1e-4;
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
    lg.Location = 'best';
    lg.Color    = 'none';
    yline(0);
    xline(0);
end
%% Salvataggio

if flag ==1
    save(1:3 ,'car','immagini/Weick_complete/caratteristiche/');
    save(4:15,'distr','immagini/Weick_complete/distribuzioni/');
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
    save(16:18 ,'car_comp','immagini/Weick_complete/caratteristiche/');
end

% % MODEL 3D
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
        k  = 6; k2 = 0.7;
        Cd = (Cd-k2) + k2*(1+0.25*M^k);
    else
        Cd = NaN;
    end
end



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