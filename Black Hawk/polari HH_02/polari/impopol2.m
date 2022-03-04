%% Importa polari profilo HH_02
clc; clear; close all
ftsize = 10;
m2tflag = 0;
d = ls('HU*');
Relist = [1250;2500;5000;10000;500];
formatspec={'-','--','-.',':'};
formatspec2={'.','*','s','^'};

for i = 1:length(Relist)-1
    nomefile       = d(i,:);
    data           = importdata(nomefile,' ',12);
    aero.alpha     = data.data(:,1);
    aero.Cl        = data.data(:,2);
    aero.Cd        = data.data(:,3);
    aero.Cl_a      = (aero.Cl(2:end)-aero.Cl(1:end-1)).*...
                    (aero.alpha(2:end)-aero.alpha(1:end-1)).^-1;
    aero.Cl_a_mode = mode(aero.Cl_a);
    save(['Aero_HH02_Re',num2str(Relist(i)),'.mat'],'aero')
    
    figure(1)
    plot(aero.alpha,aero.Cl,[formatspec{i},'k'],'DisplayName',['Re = ',num2str(Relist(i))])
    xlabel('$\alpha$','Interpreter','Latex','FontSize',ftsize);
    ylabel('$C_l$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
    hold on

    figure(2)
    plot(aero.alpha,aero.Cd,[formatspec{i},'k'],'DisplayName',['Re = ',num2str(Relist(i))])
    hold on
    xlabel('$\alpha$','Interpreter','Latex','FontSize',ftsize);
    ylabel('$Cd$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);

    figure(3)
    plot(aero.alpha(1:end-1),aero.Cl_a,[formatspec2{i},'k'],'DisplayName',['Re = ',num2str(Relist(i))])
    xlabel('$\alpha$','Interpreter','Latex','FontSize',ftsize);
    ylabel('$C_{l,\alpha} [1/^{circ}]$','Interpreter','Latex','FontSize',ftsize,'Rotation',90);
    hold on

    plot(aero.alpha(1:end-1),ones(length(aero.alpha)-1,1)*aero.Cl_a_mode,[formatspec{i},'k'],'DisplayName',['mode = ',num2str(i)])
    ylim([0 3*pi*(pi/180)])
end
for i = 1:3
    figure(i)
    grid on
    ax = gca;
    ax.FontSmoothing = 'on';
    ax.TickLabelInterpreter = 'latex';
    ax.TickLength = [0.005 0.025];
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    leg = legend();
    leg.AutoUpdate  = 'off';
    leg.NumColumns  = 1;
    leg.Location    = 'best';
    leg.Orientation = 'vertical';
    leg.Interpreter = 'latex';
    leg.Color = 'none';
end

h_fig_1 = figure(1);
h_fig_2 = figure(2);
h_fig_3 = figure(3);

if m2tflag == 1
        matlab2tikz('filename','cl_vs_alpha', 'figurehandle', h_fig_1);
        matlab2tikz('filename','cd_vs_alpha', 'figurehandle', h_fig_2);
        matlab2tikz('filename','cl_alpha', 'figurehandle', h_fig_3);
    else
        disp('Non stai generando nessun file .tex!');
end
