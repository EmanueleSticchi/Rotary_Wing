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
V_inf    = linspace(0,20);           % Climb speed [m/s]
rotore.h = 0;                        % Assume density air = 1.23 Kg/m^3
rotore   = rotore.ambient();         % compute ambient conditions
Vtheta0_d= linspace(10,19,100);
Vtheta0  = convang(Vtheta0_d,'deg','rad');                    % root pitch (comando collettivo) in [1,19 deg]

%% Analysis ---------------------------------------------------------------
options   = BEMTset_rotor();          
options.B = 0.97;                    % set B for tip correction
formatspec={'-','--','-.',':'};
ik = 0;
for i =1:length(Vtheta0)
    theta0 = Vtheta0(i);
    rotore = rotore.BEMT_salita(V_inf,theta0,options);

    
%% Graphics
    s = rotore.Analisi_salita{rotore.n_analisi_salita,1};
    
    if mod(i,25) == 0
        ik = ik + 1;
        figure(1)
        plot(s.mu,s.Tc,[formatspec{ik},'k'])
        xlabel('$\mu$','Interpreter','Latex','FontSize',ftsize);
        ylabel('$T_c$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
        hold on

        figure(2)
        plot(s.mu,s.Qc,[formatspec{ik},'k'])
        xlabel('$\mu$','Interpreter','Latex','FontSize',ftsize);
        ylabel('$Q_c$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
        hold on
    end
    %% Figure of Merit
    Tc_hover(i) = s.Tc(s.mu == 0);
    Qc_hover(i) = s.Qc(s.mu == 0);
    FM(i)       = s.Tc(s.mu == 0)^1.5/sqrt(2)/s.Qc(s.mu == 0);
    
%     disp(['Figure of Merit, FM = ',num2str(FM(i))])
end
figure

plot(Tc_hover,Qc_hover,'-k')
xlabel('$T_c$','Interpreter','Latex','FontSize',ftsize);
ylabel('$Q_c$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';



figure

plot(Tc_hover,FM,'-k')
xlabel('$Tc$','Interpreter','Latex','FontSize',ftsize);
ylabel('$FM$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

% Ricollocazione delle figure in unica figura

for i = 1:2
    f=figure(i);
%     f.Position = [511,366,995,625];
%     lg = legend('AutoUpdate','off');
%     set(lg,...
%         'Position',[0.406784353661566 0.833685915547491 0.483870960982336 0.0714285694973458],...
%         'Orientation','horizontal',...
%         'NumColumns',5,...
%         'AutoUpdate','off');
    %    xline(0)
    %    yline(0)
    grid on
    ax = gca;
    ax.FontSmoothing = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.TickLength = [0.005 0.025];
    ax.TickDir = 'in'; 
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    leg = legend(['$\theta_0 $ = ', sprintf('%0.2f',(Vtheta0_d(25))) ,'$^{\circ}$'],...
                 ['$\theta_0 $ = ', sprintf('%0.2f',(Vtheta0_d(50))) ,'$^{\circ}$'],...
                 ['$\theta_0 $ = ', sprintf('%0.2f',(Vtheta0_d(75))) ,'$^{\circ}$'],...
                 ['$\theta_0 $ = ', sprintf('%0.2f',(Vtheta0_d(end))) ,'$^{\circ}$'],...
        'Location','northeast',...
        'Orientation','vertical',...
        'NumColumns',1,...
        'AutoUpdate','off');
    leg.Orientation = 'vertical';
    leg.Interpreter = 'latex';
    leg.Color = 'none';
    ax.XAxis.Exponent = -2;
    ax.YAxis.Exponent = -2;
end
h_fig_1 = figure(1); 
h_fig_2 = figure(2); 
%% Calcoli con Prandtl correction
formatspec={':','-','--','.-'};
options   = BEMTset_rotor();          
options.P_correction = 'on';
h_fig_hover_dTc_PvsB = figure;
ik = 0;
for i = 1:length(Vtheta0)
    if mod(i,25) == 0
        ik = ik + 1;
        index(ik) = i; 
        if ik ~= 1
        theta0 = Vtheta0(i);
        rotore = rotore.BEMT_salita(V_inf,theta0,options);
        s1     = rotore.Analisi_salita{i,1};
        s      = rotore.Analisi_salita{rotore.n_analisi_salita,1};
      
        plot(s.mu, ((s.Tc-s1.Tc)./s.Tc)*100 ,[formatspec{ik},'k']);
        hold on
        end
    end
    
end


xlabel('$\mu$','Interpreter','Latex','FontSize',ftsize);
ylabel('$\Delta T_c \%$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
% ylim([-10 10]);
grid on
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
leg = legend(['$\theta_0 $ = ', sprintf('%0.1f',Vtheta0_d(index(1))) ,'$^{\circ}$'],...
            ['$\theta_0 $ = ', sprintf('%0.1f',Vtheta0_d(index(2))) ,'$^{\circ}$'],...
            ['$\theta_0 $ = ', sprintf('%0.1f',Vtheta0_d(index(3))) ,'$^{\circ}$'],...
            ['$\theta_0 $ = ', sprintf('%0.1f',Vtheta0_d(index(4))) ,'$^{\circ}$'],...             
            'Location','north');
        leg.Orientation = 'vertical';
        leg.Interpreter = 'latex';
        leg.Color = 'none';

%% Salvataggio dei file '.tex'
if m2tflag == 1
%         matlab2tikz('filename','mu_Tc', 'figurehandle', h_fig_1);
%         matlab2tikz('filename','mu_Qc', 'figurehandle', h_fig_2);
%         matlab2tikz('filename','hover_polar', 'figurehandle', h_fig_3);
%         matlab2tikz('filename','figure_of_merit', 'figurehandle', h_fig_4);
        matlab2tikz('filename','h_fig_hover_dTc_PvsB', 'figurehandle', h_fig_hover_dTc_PvsB);
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


