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
h = 0;
%% Level flight & Perfomance
V_inf_Vec = linspace(0,convvel(350,'km/h','m/s'),50);
W_vec     = linspace(W_min,heli_1.W_mtow,70);
W_kg      = W_vec./9.81;

for i = 1:length(W_vec)
W  = W_vec(i);
heli_1 = heli_1.Performance_Heli(0,2*heli_1.engine_power,W,f,heli_1.fuel_load);
sp = heli_1.PerfA{heli_1.n_PerfA,1};
V_max(i)     = sp.V_max;
V_BE(i)      = sp.V_BE;
V_BR(i)      = sp.V_BR;
Endu(i)      = sp.Endu;
Range(i)     = sp.Range;
Vd_autorot(i)= sp.Vd_autorot_min;
ROC_max(i)   = sp.ROC_max;
gamma_max(i) = sp.gamma_max;
ROD_min(i)   = sp.ROD_min;
gamma_min(i) = sp.gamma_min;
end


%% Graphics
h_fig_1 = figure;
plot(W_kg,V_max,'^k')
xlabel('$W[kg]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$V_{max}[m/s]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);

h_fig_2 = figure;
plot(W_kg,V_BE,'^k')
xlabel('$W[kg]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$V_{BE}[m/s]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);

h_fig_3 = figure;
plot(W_kg,V_BR,'^k')
xlabel('$W[kg]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$V_{BR}[m/s]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);

h_fig_4 = figure;
plot(W_kg,Endu,'-k')
xlabel('$W[kg]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$Endurance[s]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);

h_fig_5 = figure;
plot(W_kg,Range,'-k')
xlabel('$W[kg]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$Range[m]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);

h_fig_6 = figure;
plot(W_kg,Vd_autorot,'-k')
xlabel('$W[kg]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$V_{autorot}[m/s]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);

h_fig_7 = figure;
plot(W_kg,ROC_max,'-k')
xlabel('$W[kg]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$ROC_{max}$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);

h_fig_8 = figure;
plot(W_kg,convang(gamma_max,'rad','deg'),'-k')
xlabel('$W[kg]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$\gamma_{max}[^{\circ}]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);

h_fig_9 = figure;
plot(W_kg,ROD_min,'-k')
xlabel('$W[kg]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$ROD_{max}[m/s]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);

h_fig_10 = figure;
plot(W_kg,convang(gamma_min,'rad','deg'),'-k')
xlabel('$W[kg]$','Interpreter','Latex','FontSize',ftsize);
ylabel('$\gamma_{min}[^{\circ}]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);

for i =1:10
    figure(i)
    grid on
    ax = gca;
    ax.FontSmoothing = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.TickLength = [0.005 0.025];
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
end

