%% Analisi Rotore Principale HUGHES HELICOPTERS HH-02 AIRFOIL in SALITA ASSIALE
clc; clear; close all
global aero
m2tflag = 1;
ftsize = 12;
%% Data -------------------------------------------------------------------
% geometry
rotore   = Rotor();
rotore   = rotore.r(linspace(0.1,1,20));      % dominio radiale [\]
rotore.R = convlength(24,'ft','m');           % Raggio rotore   [m]
rotore.N = 4;                                 % numero di pale  [\]
rotore.c = linspace(0.53,0.53,rotore.n_r);    % corde           [m]
I_MR     = 3800;
I_MR     = convmass(I_MR,'slug','kg');
I_MR     = convlength(convlength(I_MR,'ft','m'),'ft','m');
rotore   = rotore.mass_prop('I',I_MR);        % Mom. di inerzia [Kg*m^2]
rotore.theta_t = convang(-9,'deg','rad');     % theta twist     [deg]

% aerodynamics
load('polari HH_02\polari\Aero_HH02_Re1250.mat')
rotore.Cl_alpha = aero.Cl_a_mode*180/pi;
rotore.Cl       = @(alpha) CL_(alpha);
rotore.Cd       = @(alpha) CD_(alpha);

% working conditions
rotore   = rotore.rot_vel('omega',convvel(726,'ft/s','m/s')/rotore.R);
V_inf    = linspace(0,20,4);           % Climb speed [m/s]
rotore.h = 0;                        % Assume density air = 1.23 Kg/m^3
rotore   = rotore.ambient();         % compute ambient conditions
Vtheta0  = convang([10, 13, 16, 19]...
    ,'deg','rad');                   % root pitch (comando collettivo) in [1,19 deg]

%% Analysis ---------------------------------------------------------------
options   = BEMTset_rotor();          
options.B = 0.97;                    % set B for tip correction

%% Hover & Salita assiale
for i = 1:length(Vtheta0)
    theta0 = Vtheta0(i);
    rotore = rotore.BEMT_salita(V_inf,theta0,options);
    s = rotore.Analisi_salita{rotore.n_analisi_salita,1};
end
%% Graphics
% Hover
formatspec={'-','--','.-',':'};
h_fig_hover_dTc = figure;
for i = 1:length(Vtheta0)
plot(rotore.r_bar,rotore.Analisi_salita{i,1}.dTc(1,:),[formatspec{i},'k']);
hold on
end
xlabel('$\bar{r}$','Interpreter','Latex','FontSize',ftsize);
ylabel('$dT_c/dr$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
leg = legend(['$\theta_0 $ = ', sprintf('%0.2f',(Vtheta0(1))) ,'$^{\circ}$'],['$\theta_0 $ = ', sprintf('%0.2f',(Vtheta0(2))) ,'$^{\circ}$'],...
    ['$\theta_0 $ = ', sprintf('%0.2f',(Vtheta0(3))) ,'$^{\circ}$'],['$\theta_0 $ = ', sprintf('%0.2f',(Vtheta0(4))) ,'$^{\circ}$'],'Location','north');
leg.Orientation = 'vertical';
leg.Interpreter = 'latex';
leg.Color = 'none';

h_fig_hover_dQc = figure;
for i = 1:length(Vtheta0)
plot(rotore.r_bar,rotore.Analisi_salita{i,1}.dQc(1,:),[formatspec{i},'k']);
hold on
end
xlabel('$\bar{r}$','Interpreter','Latex','FontSize',ftsize);
ylabel('$dQ_c/dr$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
leg = legend(['$\theta_0 $ = ', sprintf('%0.2f',(Vtheta0(1))) ,'$^{\circ}$'],['$\theta_0 $ = ', sprintf('%0.2f',(Vtheta0(2))) ,'$^{\circ}$'],...
    ['$\theta_0 $ = ', sprintf('%0.2f',(Vtheta0(3))) ,'$^{\circ}$'],['$\theta_0 $ = ', sprintf('%0.2f',(Vtheta0(4))) ,'$^{\circ}$'],'Location','north');
leg.Orientation = 'vertical';
leg.Interpreter = 'latex';
leg.Color = 'none';

% Salita assiale
h_fig_sas_dTc = figure;
for i = 1:length(V_inf)
plot(rotore.r_bar,rotore.Analisi_salita{3,1}.dTc(i,:),[formatspec{i},'k']);
hold on
end
xlabel('$\bar{r}$','Interpreter','Latex','FontSize',ftsize);
ylabel('$dT_c/dr$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
leg = legend(['$\mu = $',sprintf('%0.2f',(s.mu(1))),'~'],['$\mu = $',sprintf('%0.2f',(s.mu(1))),'~'],['$\mu = $',sprintf('%0.2f',(s.mu(3))),'~'],...
    ['$\mu = $',sprintf('%0.2f',(s.mu(4))),'~'],'Location','northwest');
leg.Orientation = 'vertical';
leg.Interpreter = 'latex';
leg.Color = 'none';

h_fig_sas_dQc = figure;
for i = 1:length(V_inf)
plot(rotore.r_bar,rotore.Analisi_salita{3,1}.dQc(i,:),[formatspec{i},'k']);
hold on
end
xlabel('$\bar{r}$','Interpreter','Latex','FontSize',ftsize);
ylabel('$dT_c/dr$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
leg = legend(['$\mu = $',sprintf('%0.2f',(s.mu(1))),'~'],['$\mu = $',sprintf('%0.2f',(s.mu(1))),'~'],['$\mu = $',sprintf('%0.2f',(s.mu(3))),'~'],...
    ['$\mu = $',sprintf('%0.2f',(s.mu(4))),'~'],'Location','northwest');
leg.Orientation = 'vertical';
leg.Interpreter = 'latex';
leg.Color = 'none';



%% Salvataggio dei file '.tex'
if m2tflag == 1
        matlab2tikz('filename','hover_dTc', 'figurehandle', h_fig_hover_dTc);
        matlab2tikz('filename','hover_dQc', 'figurehandle', h_fig_hover_dQc);
        matlab2tikz('filename','sas_dTc', 'figurehandle', h_fig_sas_dTc);
        matlab2tikz('filename','sas_dQc', 'figurehandle', h_fig_sas_dQc);
    else
        disp('Non stai generando nessun file .tex!');
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


