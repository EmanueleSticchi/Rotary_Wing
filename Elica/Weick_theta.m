%%  Analisi di un'elica Weick con la BEMT al varia del passo di radice
clc; clear; close all
flag = 1;
%%  DATI ----------------------------------------------------------
% importo l'oggeto elica
load el_Weick

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
