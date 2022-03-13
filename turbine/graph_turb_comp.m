%% Importa geom
clc; clear; close all
folder = 'immagini/turb/';
%
load turb
load turbina_proggettata.mat

figure; hold on
plot(turb.RNodes,turb.sChord,'k')
plot(obj.r,obj.c,'--k')
xlabel('$r [m]$','Interpreter','latex',...
    'FontSize',12)
ylabel('$c [m]$','Interpreter','latex',...
    'Rotation',90,'FontSize',12)

grid on

figure; hold on
plot(turb.RNodes,turb.AeroT,'k')
plot(obj.r,obj.theta*180/pi,'--k')
xlabel('$r [m]$','Interpreter','latex',...
    'FontSize',12)
ylabel('$\theta [deg]$','Interpreter','latex',...
    'Rotation',90,'FontSize',12)

grid on


count = 0;
for i =1:2
    count = count + 1;
    figure(i)
    FileName = sprintf(['turb_comp','%d.eps'], count);
    ax = gca;
    exportgraphics(ax,[folder,FileName])
end