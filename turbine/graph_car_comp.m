%% Turb comparison
% Confronto tra la turbina progettata con Qblade e quella progettata con il
% metodo tog (par 9.7)
clc; clear; close all
folder = 'immagini/caratteristiche/';
% importa qblade
load cp
lam = table2array(cp(:,1));
cp  = table2array(cp(:,2)); 
load ct
% importa togn
load cp_prog_togn.mat
load ct_prog_togn.mat

figure; hold on
plot(lam,cp,'k','DisplayName','Q')
cp_togn  = cpprogtogn(cpprogtogn(:,2)>0,2);
lam_togn = cpprogtogn(cpprogtogn(:,2)>0,1);
plot(lam_togn,cp_togn,'--k','DisplayName','T')
xlabel('$\lambda = \frac{\Omega R}{V_{\infty}}$',...
    'Interpreter','latex',...
    'FontSize',12)
ylabel('$C_P$','Interpreter','latex',...
    'Rotation',90,'FontSize',12)
grid on
%
lam = table2array(ct(:,1));
ct  = table2array(ct(:,2)); 

figure; hold on
plot(lam,ct,'k','DisplayName','Q')
ct_togn = ctprogtogn(cpprogtogn(:,2)>0,2);
plot(lam_togn,ct_togn,'--k','DisplayName','T')
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
%     legend
    FileName = sprintf(['c_comp','%d.eps'], count);
    ax = gca;
    exportgraphics(ax,[folder,FileName])
end