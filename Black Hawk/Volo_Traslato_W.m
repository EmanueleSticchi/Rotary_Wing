%%  Analisi Rotore Principale in volo tralato
clc; clear; close all
global aero
m2tflag = 0;
ftsize = 12;
%% Data -------------------------------------------------------------------
rotore   = Rotor();
% aerodynamics
load('polari HH_02\polari\Aero_HH02_Re1250.mat')
rotore.Cl_alpha = aero.Cl_a_mode*180/pi;
rotore.Cl       = @(alpha) CL_(alpha);
rotore.Cd       = @(alpha) CD_(alpha);

% working conditions
rotore   = rotore.rot_vel('omega',726/24); % [rad/s]
rotore.h = 0;                        % Assume density air = 1.23 Kg/m^3
rotore   = rotore.ambient();         % compute ambient conditions
V_inf    = linspace(0.1,convvel(350,'km/h','m/s'));
Chi      = convang(0.001,'deg','rad');
f        = 2;
WMTOW    = convforce(17650,'lbf','N');
vW       = [0.9, 1, 1.1];

% geometry
rotore   = rotore.r(linspace(0.1,1,20));      % dominio radiale [\]
rotore.R = convlength(24,'ft','m');           % Raggio rotore   [m]
rotore.N = 4;                                 % numero di pale  [\]
rotore.c = linspace(0.53,0.53,rotore.n_r);    % corde           [m]
I_MR     = 3800;
I_MR     = convmass(I_MR,'slug','kg');
I_MR     = convlength(convlength(I_MR,'ft','m'),'ft','m');
rotore   = rotore.mass_prop('I',I_MR);        % Mom. di inerzia [Kg*m^2]
rotore   = rotore.mass_prop('G',7);
rotore.theta_t = convang(-9,'deg','rad');     % theta twist     [deg]

formatspec = {'-k';'--k';'.-k'};
% prename = '$\frac{W}{W_{MTOW}} = ';
prename = '$W/W_{MTOW} = ';
for i =1: length(vW)
    W=vW(i)*WMTOW;
    name = [prename,num2str(vW(i)),'$'];
    %%  Analisys --------------------------------------------------------------
    rotore = rotore.BEMT_articulated('T',W,V_inf,Chi,f,BEMTset_rotor());
    
    %% Graphics ---------------------------------------------------------------
    s = rotore.Analisi_articulated{i,1};
    h_fig_mu_Pc_trasl = figure(1);
    plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.Pc_Vec,formatspec{i,1},'DisplayName',name); %!!
    xlabel('$\mu$','Interpreter','Latex','FontSize',ftsize);
    ylabel('$P_{c}$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
    grid on
    ax = gca;
    ax.FontSmoothing = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.TickLength = [0.005 0.025];
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    hold on
    
    h_fig_mu_Tc_trasl = figure(2);
    plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.Tc,formatspec{i,1},'DisplayName',name);
    xlabel('$\mu$','Interpreter','Latex','FontSize',ftsize);
    ylabel('$T_c$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
    grid on
    ax = gca;
    ax.FontSmoothing = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.TickLength = [0.005 0.025];
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    hold on
    
    h_fig_mu_theta0_trasl = figure(3);
    plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.theta0*180/pi,formatspec{i,1},'DisplayName',name);
    xlabel('$\mu$','Interpreter','Latex','FontSize',ftsize);
    ylabel('$\theta_0$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
    grid on
    ax = gca;
    ax.FontSmoothing = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.TickLength = [0.005 0.025];
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    hold on
    
    h_fig_mu_a_TPP_trasl = figure(4);
    plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.alpha_TPP_Vec*180/pi,formatspec{i,1},'DisplayName',name);
    xlabel('$\mu$','Interpreter','Latex','FontSize',ftsize);
    ylabel('$\alpha_{TPP}$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
    grid on
    ax = gca;
    ax.FontSmoothing = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.TickLength = [0.005 0.025];
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    hold on
    
    h_fig_mu_Hc_trasl = figure(5);
    plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.Hc_Vec,formatspec{i,1},'DisplayName',name);
    xlabel('$\mu$','Interpreter','Latex','FontSize',ftsize);
    ylabel('$H_c$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
    grid on
    ax = gca;
    ax.FontSmoothing = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.TickLength = [0.005 0.025];
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    hold on
    
    h_fig_mu_Yc_trasl = figure(6);
    plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.Yc_Vec,formatspec{i,1},'DisplayName',name);
    xlabel('$\mu$','Interpreter','Latex','FontSize',ftsize);
    ylabel('$Y_c$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
    grid on
    ax = gca;
    ax.FontSmoothing = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.TickLength = [0.005 0.025];
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    hold on
    
    
    h_fig_mu_lam_trasl = figure(7);
    plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.lam_Vec,'-k');
    hold on
    plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.lam_i,'--k');
%     plot(s.mu,s.lam_c,'.-k');
    xlabel('$\mu$','Interpreter','Latex','FontSize',ftsize);
    ylabel('$\lambda~\lambda_i$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
    grid on
    ax = gca;
    ax.FontSmoothing = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.TickLength = [0.005 0.025];
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    leg = legend('$\lambda$','$\lambda_i$','Location','northwest');
    leg.Orientation = 'vertical';
    leg.Interpreter = 'latex';
    leg.Color = 'none';
    hold on
    
    h_fig_mu_flapp_trasl = figure(8);
    plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.beta0_Vec*180/pi,'-k');
    hold on
    plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.beta1c_Vec*180/pi,'--k');
    plot(s.mu.*cos(s.alpha_TPP_Vec).^-1,s.beta1s_Vec*180/pi,':k');
    xlabel('$\mu$','Interpreter','Latex','FontSize',ftsize);
    ylabel('$\beta_0~\beta_{1,s}~\beta_{1,c}$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
    grid on
    ax = gca;
    ax.FontSmoothing = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.TickLength = [0.005 0.025];
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    leg = legend('$\beta_0$','$\beta_{1,c}$','$\beta_{1,s}$','Location','southwest');
    leg.Orientation = 'vertical';
    leg.Interpreter = 'latex';
    leg.Color = 'none';
end
for i =1:6
    figure(i)
    leg = legend();
    leg.Orientation = 'vertical';
    leg.Interpreter = 'latex';
    leg.Color       = 'none';
    leg.Location    = 'northwest';
    leg.FontSize    = 12;
end 
figure(7)
annotation('textarrow',[0.782291666666666 0.836979166666667],...
    [0.613736842105262 0.482105263157895],'String',{'W'});
annotation('textarrow',[0.2421875 0.281770833333333],...
    [0.201105263157895 0.329473684210526],'String',{'W'});
figure(8)
annotation('textarrow',[0.58273381294964 0.636690647482014],...
    [0.763845605700713 0.902612826603325],'String',{'W'});
annotation('textarrow',[0.739208633093526 0.705035971223025],...
    [0.619952494061758 0.5083135391924],'String',{'W'});
annotation('textarrow',[0.732014388489209 0.663669064748205],...
    [0.422802850356295 0.299287410926367],'String',{'W'});
%% Salvataggio dei file '.tex'
folder = 'Immagini\Volo_Traslato\';
if m2tflag == 1
    matlab2tikz('filename',[folder,'mu_Pc_trasl'], 'figurehandle', h_fig_mu_Pc_trasl);
    matlab2tikz('filename',[folder,'mu_Tc_trasl'], 'figurehandle', h_fig_mu_Tc_trasl);
    matlab2tikz('filename',[folder,'mu_theta0_trasl'], 'figurehandle', h_fig_mu_theta0_trasl);
    matlab2tikz('filename',[folder,'mu_a_TPP_trasl'], 'figurehandle', h_fig_mu_a_TPP_trasl);
    matlab2tikz('filename',[folder,'mu_Hc_trasl'], 'figurehandle', h_fig_mu_Hc_trasl);
    matlab2tikz('filename',[folder,'mu_Yc_trasl'], 'figurehandle', h_fig_mu_Yc_trasl);
    matlab2tikz('filename',[folder,'mu_lam_trasl'], 'figurehandle', h_fig_mu_lam_trasl);
    matlab2tikz('filename',[folder,'flapp_trasl'], 'figurehandle', h_fig_mu_flapp_trasl);
else
    disp('Non stai generando nessun file .tex!');
end


