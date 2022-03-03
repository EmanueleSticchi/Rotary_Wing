%% Importa polari profilo HH_02
clc; clear; close all
d = ls('HU*');
Relist = [500;1250;10000;2500;5000];
for i =1:length(Relist)
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
    plot(aero.alpha,aero.Cl,'DisplayName',['Re = ',num2str(Relist(i))])
    xlabel('\alpha [deg]')
    ylabel('C_l')
    hold on
    figure(2)
    plot(aero.alpha,aero.Cd,'DisplayName',['Re = ',num2str(Relist(i))])
    hold on
    xlabel('\alpha [deg]')
    ylabel('C_d')
    figure(3)
    plot(aero.alpha(1:end-1),aero.Cl_a,'.','DisplayName',['Re = ',num2str(Relist(i))])
    xlabel('\alpha [deg]')
    ylabel('C_{l_{\alpha}} [1/deg]')
    hold on
    plot(aero.alpha(1:end-1),ones(length(aero.alpha)-1,1)*aero.Cl_a_mode)
    ylim([0 3*pi*(pi/180)])
end
for i =1 :3
    figure(i)
    legend()
end