%% Importa polari profilo VR7
clc; clear; close all
flag = 1;
d = ls('BO*');
Relist = [500;1250;10000;2500;5000;7500];
ord    = [2;3;5;6];
formatspec = {'-k';'--k';'-.k';':k'};
for k =1:4
    i = ord(k);
    nomefile       = d(i,:);
    data           = importdata(nomefile,' ',12);
    aero.alpha     = data.data(:,1);
    aero.Cl        = data.data(:,2);
    aero.Cd        = data.data(:,3);
    aero.Cl_a      = (aero.Cl(2:end)-aero.Cl(1:end-1)).*...
                    (aero.alpha(2:end)-aero.alpha(1:end-1)).^-1;
    aero.Cl_a_mode = mode(aero.Cl_a);
    
    name = sprintf('Re = %0.2e',Relist(i)*1e3);
    figure(1)
    plotta(aero.alpha,aero.Cl,...
        {'$\alpha$ [deg]';'$C_l$'},...
        formatspec{k,1},name)

    figure(2)
    plotta(aero.alpha,aero.Cd,...
        {'$\alpha$ [deg]';'$C_d$'},...
        formatspec{k,1},name)
    
end
for i =1 :2
    figure(i)
    lg = legend();
    lg.Location = 'northwest';
    lg.Color = 'none';
end

if flag ==1
    save(1:2,'pol','immagini/')
end
%% Function
function save(idxF,prename,folder)
count = 0;
for i =1:length(idxF)
    count = count + 1;
    figure(idxF(i))
    %     legend
    FileName = sprintf([prename,'%d.eps'], count);
    ax = gca;
    exportgraphics(ax,[folder,FileName])
end
end