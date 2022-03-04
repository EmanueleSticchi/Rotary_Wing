clear all; close all; clc;
ftsize= 12;
m2tflag = 0;
% Helicopter: class definition
heli_1            = Helicopter();

% Rotors: class definition and rotors properties
% Main Rotor 
heli_1.MR         = Rotor();
heli_1.MR         = heli_1.MR.r(0.1:0.01:0.9);
heli_1.MR.R       = convlength(24,'ft','m'); 
heli_1.MR.N       = 4;
heli_1.MR.c       = linspace(0.53,0.53,heli_1.MR.n_r);
heli_1.MR.omega   = convvel(726,'ft/s','m/s')/heli_1.MR.R;
heli_1.MR.theta_t = convang(-9,'deg','rad');
I_MR              = 3800;
I_MR              = convmass(I_MR,'slug','kg');
I_MR              = convlength(convlength(I_MR,'ft','m'),'ft','m');
heli_1.MR         = heli_1.MR.mass_prop('I',I_MR);
% Main rotor airfoil: HH-02

% Tail Rotor 
heli_1.TR         = Rotor();
heli_1.TR         = heli_1.TR.r(0.1:0.01:0.9);
heli_1.TR.R       = convlength(4.6,'ft','m');
heli_1.TR.N       = 4;
heli_1.TR.c       = linspace(convlength(0.83,'ft','m'),convlength(0.83,'ft','m'),heli_1.TR.n_r);
heli_1.TR.omega   = convvel(677,'ft/s','m/s')/heli_1.TR.R;
heli_1.TR.theta_t = convang(-8,'deg','rad');
I_TR              = 10;
I_TR              = convmass(I_TR,'slug','kg');
I_TR              = convlength(convlength(I_TR,'ft','m'),'ft','m');
heli_1.TR         = heli_1.TR.mass_prop('I',I_TR);
heli_1.lr         = 17.73 - heli_1.TR.R - heli_1.MR.R;
% Tail rotor airfoil: NACA 63-414

% Mass
W_empty          = 5352*9.81;
heli_1.W_fuel    = 1108*9.81;
heli_1.fuel_load = 1;
heli_1.W_mtow    = convforce(21000,'lbf','N');
W_min            = W_empty + heli_1.fuel_load*heli_1.W_fuel;


% Power
heli_1.engine_power  = 1279000;
heli_1.engine_number = 2;
P_available          = heli_1.engine_number*heli_1.engine_power;
heli_1.SFC           = 0.4;
heli_1.P_req_AUX     = 20000;
heli_1.eta_t         = 1.03;
P_av = 20000;

f = 0.007;
f = f*pi*(heli_1.MR.R)^2;

%% Level flight & Perfomance
service_ceiling = 6400;
h = linspace(0,service_ceiling,4);
W = linspace(W_min,heli_1.W_mtow,4);
W_kg = W./9.81;
V_inf_Vec = linspace(0,150,50);

for i = 1:length(W)
    heli_1 = heli_1.Req_power_level_flight(0,V_inf_Vec,W(i),f,'P',2*heli_1.engine_power);
    heli_1 = heli_1.Performance_Heli(0,2*heli_1.engine_power,W(i),f,heli_1.fuel_load);
end
 

%% Graphics
formatspec={'-','--','.-',':'};

%% Massimo peso al decollo: decomposizione delle potenze e prestazioni 
h_fig_req_pow_MR = figure; 
plot(V_inf_Vec,heli_1.PA{length(W),1}.Pi_MR,'-k');
hold on;
plot(V_inf_Vec,heli_1.PA{length(W),1}.P0_MR,'--k');
plot(V_inf_Vec,heli_1.PA{length(W),1}.P_fus_MR,'.-k');
plot(V_inf_Vec,heli_1.PA{length(W),1}.P_req_MR,':k');
xlabel('$V_{\infty}[m/s]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$P[W]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
title('\textup{Rotore principale}','Interpreter','Latex','FontSize',ftsize,'Rotation',0);
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
leg = legend('$P_{i}$',...
             '$P_{0}$',...
             '$P_{\textup{fus}}$',...
             '$P_{\textup{req}}$',...
             'Location','northeast');
leg.Orientation = 'vertical';
leg.Interpreter = 'latex';
leg.Color = 'none';

h_fig_req_pow_TR = figure; 
plot(V_inf_Vec,heli_1.PA{length(W),1}.Pi_TR,'-k');
hold on;
plot(V_inf_Vec,heli_1.PA{length(W),1}.P0_TR,'--k');
plot(V_inf_Vec,heli_1.PA{length(W),1}.P_req_MR,'.-k');
xlabel('$V_{\infty}[m/s]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$P[W]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
title('\textup{Rotore anti-coppia di coda }','Interpreter','Latex','FontSize',ftsize,'Rotation',0);
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
leg = legend('$P_{i}$',...
             '$P_{0}$',...
             '$P_{\textup{req}}$',...
             'Location','northeast');
leg.Orientation = 'vertical';
leg.Interpreter = 'latex';
leg.Color = 'none';

h_fig_req_pow_tot = figure; 
plot(V_inf_Vec,heli_1.PA{length(W),1}.P_req_TR,'--k');
hold on;
plot(V_inf_Vec,heli_1.PA{length(W),1}.P_req_MR,'-k');
plot(V_inf_Vec,heli_1.P_req_AUX*ones(length(V_inf_Vec),1),'.-k');
plot(V_inf_Vec,heli_1.PA{length(W),1}.P_av*ones(length(V_inf_Vec),1),':k');
plot(V_inf_Vec,heli_1.PA{length(W),1}.P_req_hori,'LineWidth',2,'Color','k');
xlabel('$V_{\infty}[m/s]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$P[W]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
title('\textup{Potenze totali}','Interpreter','Latex','FontSize',ftsize,'Rotation',0);
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
leg = legend('$P_{req,TR}$',...
             '$P_{req,MR}$',...
             '$P_{aux}$',...
             '$P_{av}$',...
             '$P_{tot}$',...
             'Location','northeast');
leg.Orientation = 'vertical';
leg.Interpreter = 'latex';
leg.Color = 'none';


%% Variazione di peso, sola potenza richiesta totale & la velocità di salita
line_type = (linspace(1,2,4));
h_fig_req_W = figure; 
for i = 1:length(W)
plot(V_inf_Vec,heli_1.PA{i,1}.P_req_hori,[formatspec{i},'k']);
hold on
end
plot(V_inf_Vec,P_available*ones(length(V_inf_Vec),1),'Color','k','LineWidth',1.5);
for i = 1:length(W)
plot(heli_1.PerfA{i,1}.V_BE,heli_1.PerfA{i,1}.P_min,'o','MarkerSize',10,'Color','k')
hold on
plot(heli_1.PerfA{i,1}.V_BR,heli_1.PerfA{i,1}.P_BR,'s','MarkerSize',10,'Color','k')
plot(heli_1.PerfA{i,1}.V_max,P_available,'^','MarkerSize',10,'Color','k')
end
% plot(heli_1.PerfA{1,1}.V_BE,heli_1.PerfA{1,1}.P_min,'*','DisplayName','P_{min}','MarkerSize',10)
% plot(heli_1.PerfA{1,1}.V_BR,heli_1.PerfA{1,1}.P_BR,'*','DisplayName','P_{BR}','MarkerSize',10)
% plot(heli_1.PerfA{1,1}.V_max,heli_1.PerfA{1,1}.P_av,'*','DisplayName','V_{max}','MarkerSize',10)
xlabel('$V_{\infty}[m/s]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$P[W]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
title('\textup{Potenze totali}','Interpreter','Latex','FontSize',ftsize,'Rotation',0);
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
leg = legend(['$P_{tot}(W = $', sprintf('%0.0f',(W_kg(1))) '$~kg)$'],...
             ['$P_{tot}(W = $', sprintf('%0.0f',(W_kg(2))) '$~kg)$'],...
             ['$P_{tot}(W = $', sprintf('%0.0f',(W_kg(3))) '$~kg)$'],...
             ['$P_{tot}(W = $', sprintf('%0.0f',(W_kg(4))) '$~kg)$'],...
             '$p_{av}$',...
             '$V_{BE}$',...
             '$V_{BR}$',...
             '$V_{max}$',...
             'Location','southeast');
leg.Orientation = 'vertical';
leg.Interpreter = 'latex';
leg.Color = 'none';

max_ROC = max(heli_1.PA{i,1}.Vc);

h_fig_ROC = figure; 
for i = 1:length(W)
plot(V_inf_Vec,heli_1.PA{i,1}.Vc,[formatspec{i},'k']);
hold on
end
plot(V_inf_Vec, zeros(length(V_inf_Vec),1),'-k')
plot(V_inf_Vec, zeros(length(V_inf_Vec),1),'-k')

xlabel('$V_{\infty}[m/s]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$ROC[m/s]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
leg = legend(['$P_{tot}(W = $', sprintf('%0.0f',(W_kg(1))) '$~kg)$'],...
             ['$P_{tot}(W = $', sprintf('%0.0f',(W_kg(2))) '$~kg)$'],...
             ['$P_{tot}(W = $', sprintf('%0.0f',(W_kg(3))) '$~kg)$'],...
             ['$P_{tot}(W = $', sprintf('%0.0f',(W_kg(4))) '$~kg)$'],...
             'Location','south');
leg.Orientation = 'vertical';
leg.Interpreter = 'latex';
leg.Color = 'none';



