%% Costruzione di look-up table per i dati aerodinamici
clc; clear; close all
filename = ls('Ae*');
Relist = [10000;1250;2500;500;5000]*1e3;
for i =1 : size(filename)
    load(filename(i,:))
    aerod{i,1} = aero; 
    aerod{i,2} = Relist(i);
end
v_alpha = linspace(-20,20,400);   
p = 0.999999;                           % smoothing parameter
for i=1:length(aerod)
    pp_Cl_alpha = csaps(aerod{i,1}.alpha, aerod{i,1}.Cl, p); % smoothing spline
    m_Cl(:,i) = fnval(pp_Cl_alpha, v_alpha); % funzione interpolante
end
[m_alpha,m_Re] = meshgrid(v_alpha,Relist);



interp2(m_alpha,m_Re,m_Cl,10,5e5)


