%% Importa profili Elica Weick
clc; clear; close all
c    = [6.97,...     % Corda delle sezioni, [in]
    7.36,8.08,8.41,8.43,8.17,7.53,6.53,5.21,3.74];
data = importdata('Airfoil_Weick.txt');
AF.S = data(1:11)/100;
AF.c = c;
n    = 10;
k    = 12;
for i = 1:3
    AF.up(:,i)  = data(k:k+10);
    k = k + 11;
    AF.low(:,i) = - data(k:k+10);
    k = k + 11;
    
    figure
    plot(AF.S*c(i),AF.up(:,i),'k',AF.S*c(i),AF.low(:,i))
    daspect([1 1 1])
    grid on 
    xlabel('x [in]')
    ylabel('z [in]')
end

for i = 4:n
    AF.up(:,i)  = data(k:k+10);
    k = k + 11;
    
    figure
    plot(AF.S*c(i),AF.up(:,i),'k')
    hold on
    yline(0,'k')
    daspect([1 1 1])
    grid on 
    xlabel('x [in]')
    ylabel('z [in]')
end

%% Plot dei profili normalizzati
figure
hold on
for i = 1:3    
    plot(AF.S,AF.up(:,i)/c(i),AF.S,AF.low(:,i)/c(i))
end
daspect([1 1 1])
grid on 
xlabel('x/c')
ylabel('z/c')

figure
hold on
yline(0,'k')
for i = 4:n
    plot(AF.S,AF.up(:,i)/c(i))
end
daspect([1 1 1])
grid on 
xlabel('x/c')
ylabel('z/c')