%% Analisi Rotore Principale HUGHES HELICOPTERS HH-02 AIRFOIL in SALITA ASSIALE
clc; clear; close all
global aero
%% Data -------------------------------------------------------------------
% geometry
rotore   = Rotor();
rotore   = rotore.r(linspace(0.1,1,20));      % dominio radiale [\]
rotore.R = 14.63;                             % Raggio rotore   [m]
rotore.N = 4;                                 % numero di pale  [\]
rotore.c = linspace(0.53,0.53,rotore.n_r);    % corde           [m]
I_MR     = 3800;
I_MR     = convmass(I_MR,'slug','kg');
I_MR     = convlength(convlength(I_MR,'ft','m'),'ft','m');
rotore   = rotore.mass_prop('I',I_MR);        % Mom. di inerzia [Kg*m^2]
rotore.theta_t = convang(-9,'deg','rad');     % theta twist     [deg]

% aerodynamics
load('polari HH_02\polari\Aero_HH02_Re1250.mat')
rotore.Cl_alpha = aero.Cl_a_mode;
rotore.Cl       = @(alpha) CL_(alpha);
rotore.Cd       = @(alpha) CD_(alpha);

% working conditions
rotore   = rotore.rot_vel('omega',convvel(726,'ft/s','m/s')/rotore.R);
V_inf    = linspace(0,10);           % Climb speed [m/s]
rotore.h = 0;                        % Assume density air = 1.23 Kg/m^3
rotore   = rotore.ambient();         % compute ambient conditions
Vtheta0  = convang([10,12,15]...
    ,'deg','rad');                   % root pitch (comando collettivo) in [1,19 deg]

%% Analysis ---------------------------------------------------------------
options   = BEMTset_rotor();          
options.B = 0.97;                    % set B for tip correction
for i =1:length(Vtheta0)
    theta0 = Vtheta0(i);
    rotore = rotore.BEMT_salita(V_inf,theta0,options);

%% Graphics
    s = rotore.Analisi_salita{rotore.n_analisi_salita,1};
    name = ['\theta_0 = ',num2str(theta0*180/pi),' deg'];
    
    figure(1)
    plot(s.mu,s.Tc,'DisplayName',name)
    xlabel('\mu')
    ylabel('T_c')
    hold on 

    figure(2)
    plot(s.mu,s.Qc,'DisplayName',name)
    xlabel('\mu')
    ylabel('Q_c')
    hold on 

    figure(3)
    plot(s.Tc,s.Qc,'DisplayName',name)
    xlabel('T_c')
    ylabel('Q_c')
    hold on 

    %% Figure of Merit

    FM(i) = s.Tc(s.mu == 0)^1.5/sqrt(2)/s.Qc(s.mu == 0);
    disp(['Figure of Merit, FM = ',num2str(FM(i))])
end
figure
plot(Vtheta0*180/pi,FM)
xlabel('\theta_0 [deg]')
ylabel('FM')
for i =1 : 3
   figure(i) 
   legend('AutoUpdate','off')
   xline(0)
   yline(0)
   grid on
end
%% Function ---------------------------------------------------------------
function CL=CL_(alpha)
global aero

    alpha=alpha*180/pi;
    if alpha <  aero.alpha(1)
        CL=aero.Cl(1);
    elseif alpha > aero.alpha(end)
        CL=aero.Cl(end);
    else 
        CL=interp1(aero.alpha,aero.Cl,alpha);
    end
end

function CD=CD_(alpha)
global aero
    alpha=alpha*180/pi;
    if alpha <  aero.alpha(1)
        CD=aero.Cd(1);
    elseif alpha > aero.alpha(end)
        CD=aero.Cd(end);
    else 
        CD=interp1(aero.alpha,aero.Cd,alpha);
    end
end