% Rappresentazione profilo
clear all; close all; clc;
ftsize = 10;
m2tflag = 1;

airfoil_input_table = readtable('HH-02');
airfoil_input = table2array(airfoil_input_table);
dim = length(airfoil_input(:,1));
airfoil_chord_u = airfoil_input(1:(dim/2), 1);
airfoil_chord_l = airfoil_input(((dim/2)+1):dim, 1);
airfoil_upper = airfoil_input(1:((dim)/2), 2);
airfoil_lower = airfoil_input(((dim/2)+1):dim, 2);
airfoil_lower(1) = airfoil_upper(end);
airfoil_meanline = 0.5*(airfoil_upper + airfoil_lower);

%% Post-processing e grafica
%-------------------------------------------------------------------------%
% Rappresentazione del profilo
%-------------------------------------------------------------------------%
h_fig_profilo = figure;
plot(airfoil_chord_u(:,1), airfoil_upper(:,1), '-k');
hold on;
plot(airfoil_chord_l(:,1), airfoil_lower(:,1), '-k');
title('\textbf{HH-02}','Interpreter','Latex','FontSize',ftsize);
xlabel('$x/c$','Interpreter','Latex','FontSize',ftsize);
ylabel('$z/c$','Interpreter','Latex','FontSize',ftsize,'Rotation',0,'Margin',10,'Position',[-0.07 -0.0015 0]);
grid on;
axis equal
axis([-0.005 1.005 -0.1 0.12]);
ax = gca;
ax.FontSmoothing = 'on';
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0.005 0.025];
ax.TickDir = 'in';
ax.YTick = [-0.10 0 0.10];

if m2tflag == 1
        matlab2tikz('filename','airfoil_MR', 'figurehandle', h_fig_profilo);
    else
        disp('Non stai generando nessun file .tex!');
end