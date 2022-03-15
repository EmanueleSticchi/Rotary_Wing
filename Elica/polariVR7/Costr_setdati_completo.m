%% Costruzione di look-up table per i dati aerodinamici
clc; clear; close all
global m_alpha m_Re m_Cl m_Cd

filename = ls('Ae*Re*');
Relist   = [10000;1250;2500;500;5000;7500]*1e3;
Re       = sort(Relist);
ord      = [6;2;3;1;4;5];
aerod    = cell(length(Re),2);
for i =1 : size(filename)
    load(filename(i,:))
    if aero.alpha(1) > -30
        aero.alpha = [-30;aero.alpha];
        aero.Cl    = [aero.Cl(1);aero.Cl];
        aero.Cd    = [aero.Cd(1);aero.Cd];
    end
    if aero.alpha(end) < 30
        aero.alpha = [aero.alpha;30];
        aero.Cl    = [aero.Cl;aero.Cl(end)];
        aero.Cd    = [aero.Cd;aero.Cd(end)];
    end
    aerod{ord(i),1} = aero; 
    aerod{ord(i),2} = Relist(i);
end
v_alpha = linspace(-30,30,400);   
p = 0.999999;                           % smoothing parameter
for i=1:length(aerod)
    pp_Cl_alpha = csaps(aerod{i,1}.alpha, aerod{i,1}.Cl, p); % smoothing spline
    pp_Cd_alpha = csaps(aerod{i,1}.alpha, aerod{i,1}.Cd, p); % smoothing spline
    m_Cl(i,:) = fnval(pp_Cl_alpha, v_alpha); % funzione interpolante
    m_Cd(i,:) = fnval(pp_Cd_alpha, v_alpha); % funzione interpolante
end
[m_alpha,m_Re] = meshgrid(v_alpha,Re);

Cl_(pi/180*9.95,1.2,100e5)
Cd_(pi/180*9.95,0.72,100e5)
M = linspace(0,1,200).^0.5;
for i=1:length(M) 
    Cd(i) = Cd_(pi/180*9.95,M(i),100e5); 
end
plot(M,Cd)
axis equal

% save('AeroVR7_complete.mat','m_alpha','m_Re','m_Cl','m_Cd')


%% Function

function Cl = Cl_(alpha,M,Re)
% angoli in radianti
    global m_alpha m_Re m_Cl
    alpha     = convang(alpha,'rad','deg');
    alpha_lim = [m_alpha(1,1),m_alpha(1,end)];
    Re_lim    = [m_Re(1,1)   ,m_Re(end,1)];
    
    Cl = interp2(m_alpha,m_Re,m_Cl,...
        min(max(alpha_lim(1),alpha),alpha_lim(2)),...
        min(max(Re_lim(1),Re),Re_lim(2)));
    if M <1
        Cl = Cl.*sqrt(1-M.^2).^-1;
    else
        Cl = NaN;
    end
end

function Cd = Cd_(alpha,M,Re)
% angoli in radianti
    global m_alpha m_Re m_Cd
    alpha     = convang(alpha,'rad','deg');
    alpha_lim = [m_alpha(1,1),m_alpha(1,end)];
    Re_lim    = [m_Re(1,1)   ,m_Re(end,1)];
    
    Cd = interp2(m_alpha,m_Re,m_Cd,...
        min(max(alpha_lim(1),alpha),alpha_lim(2)),...
        min(max(Re_lim(1),Re),Re_lim(2)));
    if M <1
        k  = 6; k2 = 0.7;
        Cd = (Cd-k2) + k2*(1+0.25*M^k);
    else
        Cd = NaN;
    end
end

