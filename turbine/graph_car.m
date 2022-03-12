%% Importa geom
clc; clear; close all
folder = 'immagini/caratteristiche/';
%
load cp
lam = table2array(cp(:,1));
cp  = table2array(cp(:,2)); 

figure
plot(lam(:),cp(:),'k')
xlabel('$\lambda = \frac{\Omega R}{V_{\infty}}$',...
    'Interpreter','latex',...
    'FontSize',12)
ylabel('$C_P$','Interpreter','latex',...
    'Rotation',90,'FontSize',12)
grid on
%
load ct
lam = table2array(ct(:,1));
ct  = table2array(ct(:,2)); 

figure
plot(lam,ct,'k')
xlabel('$\lambda = \frac{\Omega R}{V_{\infty}}$',...
    'Interpreter','latex',...
    'FontSize',12)
ylabel('$C_T$','Interpreter','latex',...
    'Rotation',90,'FontSize',12)
grid on


count = 0;
for i =1:2
    count = count + 1;
    figure(i)
    FileName = sprintf(['c','%d.eps'], count);
    ax = gca;
    exportgraphics(ax,[folder,FileName])
end
