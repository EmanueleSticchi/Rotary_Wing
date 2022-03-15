%% Progetto Elica Lenta per il velivolo Tecnam P92
% assumiamo profilo alare Boeing Vertol VR7
clc; clear; close all
flag = 0;
global m_alpha m_Re m_Cl m_Cd

%% DATI -------------------------------------------------------------------
% velivolo
S   = 13.2;      % m^2
b   = 8.7;
AR  = b^2/S;
e   = 0.85;      % assunto da noi
W   = 450*9.81;  % N
Vc  = 51.43;     % m/s da Janes
CD0 = 0.035;     % uguale a CD0 auunto da noi 
% elica
el    = Elica();
el.R  = 1.66/2;     % m, desunta dal trittico
r_hub = el.R*0.1; % assunto da noi
el    = el.r_(linspace(r_hub/el.R,1));
el.N  = 2;        % scelto da noi
el    = el.derived_properties();
% Aerodinamici
load AeroVR7_complete.mat
el.Cl=@(alpha,r_bar,M,Re) Cl_(alpha,M,Re);
el.Cd=@(alpha,r_bar,M,Re) Cd_(alpha,M,Re);

% funzionamento
h  = 0;
el = el.altitude(h);
el = el.rot_vel('omega',0.7*el.sound_vel/el.R);


%% scelta dell' alpha_id e Cl_id
[r,c] = size(m_alpha);

for i =1:r
%     figure(1)
%     plot(m_alpha(i,:),m_Cl(i,:))
%     hold on
%     figure(2)
%     plot(m_alpha(i,:),m_Cd(i,:))
%     hold on
%     figure(3)
%     plot(m_alpha(i,:),m_Cl(i,:).*m_Cd(i,:).^-1)
%     hold on
    [Emax,id_max] = max(m_Cl(i,:).*m_Cd(i,:).^-1);
    E_max(i) =  Emax;
    a_id(i)  = m_alpha(i,id_max);
end

% figure 
% plot(m_Re(E_max < 1e4,1),E_max(E_max < 1e4))
% figure
% plot(m_Re(E_max < 1e4,1),a_id(E_max < 1e4))

alpha_id = mean(a_id(E_max < 1e4));
for i =1:r
   [a,id] = min(abs(m_alpha(i,:)-alpha_id));
   Cl_id(i) = m_Cl(i,id);
end
Cl_id    = mean(Cl_id)*ones(el.n_r,1);
% Cl_id    = 0.2*ones(el.n_r,1);
alpha_id = alpha_id*ones(el.n_r,1)*pi/180;     

%% Progetto alla velocità di crociera
CL = W/(0.5*el.rho*Vc^2*S);     
CD = CD0 +CL^2/(pi*AR*e);
T  = 0.5*el.rho*Vc^2*S*CD;
% CT = T/(el.rho*el.n^2*(2*el.R)^4);
J_des  = Vc/(el.n*2*el.R);     
options=BEMTset();  
options.P_correction   = 'on';
options.Hub_correction = 'on';
options.Cd_hub         = 0.5;
options.Design         = 'on';
options.Freccia_opt    = 'off';

el = el.Design(J_des,T,Cl_id,alpha_id,options);



%% Curve Caratteristiche

alpha0=-2*pi/180;
options=BEMTset();              options.P_correction = 'on';
options.Hub_correction = 'on';  options.Cd_hub       = 0.5;
options.toll = 1e-6;
v_J = 0.15:0.05:1.1;
el=el.BEMT(v_J,alpha0,options);

s = el.Analisi{1, 1};
figure(1)
plotta(s.J,s.CT,{'J = $ \frac{V_{\infty}}{nD}$';'$C_T$ = $\frac{T}{\rho n^2 D^4}$'})

figure(2)
plotta(s.J,s.CP,{'J = $ \frac{V_{\infty}}{nD}$';'$C_P$ = $\frac{P}{\rho n^3 D^5}$'})

figure(3)
log=s.CT>=0;
plotta(s.J(log),s.eta(log,1),{'J = $ \frac{V_{\infty}}{nD}$';'$\eta$ = $\frac{TV_{\infty}}{P}$'})

%% Distribuzioni in crociera
el=el.BEMT(J_des,alpha0,options);

s = el.Analisi{2, 1};

figure(4)
plotta(el.r_bar,s.dCt_dr_bar,...
    {'$ \bar{r}$';'$\frac{dC_T}{d\bar{r}}$'})

figure(5)
plotta(el.r_bar,s.dCq_dr_bar,...
    {'$ \bar{r}$';'$\frac{dC_Q}{d\bar{r}}$'})

figure(6)
plotta(el.r_bar,s.dCp_dr_bar,...
    {'$ \bar{r}$';'$\frac{dC_P}{d\bar{r}}$'})

figure(7)
plotta(el.r_bar,s.eta_e,{'$ \bar{r}$';'$\eta_e$'})


figure(8)
plotta(el.r_bar,s.Mach,{'$ \bar{r}$';'$M$'})

figure(9)
plotta(el.r_bar,s.Re,{'$ \bar{r}$';'$Re$'})

figure(10)
plotta(el.r_bar,s.alpha*180/pi,{'$ \bar{r}$';'$\alpha$'})

figure(11)
plotta(el.r_bar,s.phi*180/pi,{'$ \bar{r}$';'$\varphi$'})

figure(12)
plotta(el.r_bar,s.a,{'$ \bar{r}$';'$a$'})

figure(13)
plotta(el.r_bar,s.ap,{'$ \bar{r}$';'$a$'''})

figure(14)
plotta(el.r_bar,s.lambda1,{'$ \bar{r}$';'$\lambda_1$'})

figure(15)
plotta(el.r_bar,s.lambda2,{'$ \bar{r}$';'$\lambda_2$'})

%% Salvataggio

if flag ==1
    save(1:3 ,'car','immagini/Design/p92/caratteristiche/');
    save(4:15,'distr','immagini/Design/p92/distribuzioni/');
end

%% Graph geom (see Geom_p92.m)
data = importdata('BOEING_VERTOL_VR-7.dat');
x    = data.data(:,1);   z = data.data(:,2);
el.Model3D(x,z)
axis off
figure
plotta(el.r_bar,el.c,{'$ \bar{r}$';'c [m]'});
figure
plotta(el.r_bar,(el.theta)*180/pi,{'$ \bar{r}$';'$\theta$ [deg]'});
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
