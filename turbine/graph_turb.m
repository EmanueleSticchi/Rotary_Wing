%% Importa geom
clc; clear; close all
folder = 'immagini/turb/';
%
load turb

figure
plot(turb.RNodes,turb.sChord,'k')
xlabel('$r [m]$','Interpreter','latex',...
    'FontSize',12)
ylabel('$c [m]$','Interpreter','latex',...
    'Rotation',90,'FontSize',12)

grid on

figure
plot(turb.RNodes,turb.AeroT,'k')
xlabel('$r [m]$','Interpreter','latex',...
    'FontSize',12)
ylabel('$\theta [deg]$','Interpreter','latex',...
    'Rotation',90,'FontSize',12)

grid on


count = 0;
for i =1:2
    count = count + 1;
    figure(i)
    FileName = sprintf(['turb','%d.eps'], count);
    ax = gca;
    exportgraphics(ax,[folder,FileName])
end