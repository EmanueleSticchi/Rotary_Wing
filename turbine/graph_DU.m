%% Importa airfoil
clc; clear; close all
folder = 'immagini/DU';
%
data = importdata('DU84-132V3.dat');
x    = data.data(:,1);
z    = data.data(:,2);

figure
plot(x,z,'k')
xlabel('$x/c$','Interpreter','latex',...
    'FontSize',12)
ylabel('$y/c$','Interpreter','latex',...
    'Rotation',90,'FontSize',12)
daspect([1 1 1])
grid on
%% Importa polari

data  = importdata('DU84-132V3_T1_Re1.000_M0.00_N9.0_360_M.dat',' ',14);
alpha = data.data(:,1);
Cl    = data.data(:,2);
Cd    = data.data(:,3);

q = {Cl;Cd;Cl.*Cd.^-1};
name = {'$C_l$';'$C_d$';'$\frac{C_l}{C_d}$'};
for i =1 :3
    figure
    plot(alpha,q{i,1},'k')
    hold on
    plot(alpha(168:856),q{i,1}(168:856),'.k')
    xlabel('$\alpha$ [deg]',...
        'Interpreter','latex','FontSize',12)
    ylabel(name{i,1},...
        'Interpreter','latex','Rotation',90,...
        'FontSize',12)
    
    grid on
end
count = 0;
for i =1:4
    count = count + 1;
    figure(i)
    FileName = sprintf(['DU','%d.eps'], count);
    ax = gca;
    exportgraphics(ax,[folder,FileName])
end

