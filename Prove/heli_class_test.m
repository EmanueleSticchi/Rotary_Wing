clear all; close all; clc;
% Helicopter: class definition
heli_1            = Helicopter();

% Rotors: class definition and rotors properties
% Main Rotor 
heli_1.MR         = Rotor();
heli_1.MR         = heli_1.MR.r(0.1:0.01:0.9);
heli_1.MR.R       = 14.63;
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
heli_1.TR.R       = 2.79;
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
heli_1.W_mtow    = W_empty + heli_1.fuel_load*heli_1.W_fuel;


% Power
heli_1.engine_power  = 1279000;
heli_1.engine_number = 2;
P_available          = heli_1.engine_number*heli_1.engine_power;
heli_1.SFC           = 0.4;
heli_1.P_req_AUX     = 20000;
heli_1.eta_t         = 1.2;

f = 0.007;
f = f*pi*(heli_1.MR.R)^2;
h = 0;
V_inf_Vec = linspace(0,convvel(200,'km/h','m/s'),50);
heli_1 = heli_1.Req_power_level_flight(h,V_inf_Vec,heli_1.W_mtow,f,'P',20000);

%% Graphics
s = heli_1.PA{1,1};
plot(V_inf_Vec,s.Vc);
