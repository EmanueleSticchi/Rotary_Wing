%%  Analisi dell'elica progettata per il velivolo c27j
% al variare del passo theta0
% modello aerodinamico completo profilo BV VR7
clc; clear; close all
flag = 0;
global m_alpha m_Re m_Cl m_Cd

%%  DATI ----------------------------------------------------------
% importo l'oggeto elica
load el_c27j.mat
% Aerodinamici
load AeroVR7_complete.mat
el.Cl=@(alpha,r_bar,M,Re) Cl_(alpha,M,Re);
el.Cd=@(alpha,r_bar,M,Re) Cd_(alpha,M,Re);

% Di funzionamento
vtheta = convang(linspace(10,50,5),'deg','rad');
el=el.altitude(0);
options=BEMTset();             options.P_correction='on'; 
options.Hub_correction = 'on'; options.Cd_hub =0.5;
% alpha0 = [-2 -10 -10 -30 -50]*pi/180;
alpha0 = 5*pi/180;
J_end  = [3 3 3 3 3];
for tdx = 1:length(vtheta)
    theta_75 = vtheta(tdx);
    [~,idx] = min(abs(el.r_bar - 0.75));
    el.theta = el.theta - el.theta(idx) + theta_75;  % setto il calettamento nominale al 75 %
    J=[0.5:0.05:J_end(tdx)];
    %% Analisi BEMT --------------------------------------------------------       
    el=el.BEMT(J,alpha0,options);
    
end

%% Post - Processing    
formatspec = {'-k';'--k';':k';'.-k';'^-k';'s-k'};
for tdx = 1:length(vtheta)
    s = el.Analisi{tdx, 1};
    name = ['\theta_{75} = ',num2str(vtheta(tdx)*180/pi),'°'];
    figure(1)
    plotta(s.J,s.CT,...
        {'J = $ \frac{V_{\infty}}{nD}$';'$C_T$ = $\frac{T}{\rho n^2 D^4}$'},...
        formatspec{tdx,1},name)
    
    figure(2)
    plotta(s.J,s.CP,...
        {'J = $ \frac{V_{\infty}}{nD}$';'$C_P$ = $\frac{P}{\rho n^3 D^5}$'},...
        formatspec{tdx,1},name)
    
    figure(3)
    log=s.CT>=0;
    plotta(s.J(log),s.eta(log,1),...
        {'J = $ \frac{V_{\infty}}{nD}$';'$\eta$ = $\frac{TV_{\infty}}{P}$'},...
        formatspec{tdx,1},name)    
end

for i =1:3
    figure(i)
    lg = legend();
    lg.AutoUpdate='off';
    lg.Location = 'best';
    lg.Color    = 'none';
    yline(0);
    xline(0);
end

if flag ==1
    save(1:3,'car_theta','immagini/Design/c27j/caratteristiche/')
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