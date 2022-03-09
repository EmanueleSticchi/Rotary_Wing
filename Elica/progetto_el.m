%%  Analisi di un'elica Weick con la BEMT
clc; clear; close all
global aero
%%  DATI 
% numerici e di controllo
toll=1e-6;
el=Elica();
options=BEMTset();      
options.Freccia_opt='off';
% geometrici
el=el.r_((0.1:0.01:1)');
el.D = 1.66; % Raggio dell'elica, [m]
el.N = 4;    % Numero di pale    
el.R=0.5*el.D;
% di funzionamento
el=el.rot_vel('RPM',2000); 
el=el.altitude(0);
% aerodinamici
% el.Cl = @(alpha,r_bar,M,Re) 2*pi*alpha;
% el.Cd = @(alpha,r_bar,M,Re) 0.02*alpha./alpha;

load Aero_NACA16212.mat
el.Cl=@(alpha,r_bar,M,Re) CL_(alpha);
el.Cd=@(alpha,r_bar,M,Re) CD_(alpha);

% if isequal(options.Freccia_opt,'off')
% %     el.LAMBDA = linspace(0,20*pi/180,el.n_r);
%     el.LAMBDA = zeros(el.n_r,1);
%     for i = 80:91
%     el.LAMBDA(i) = 0.3;
%     end
% end

%% Assegno i valori per il progetto e risolvo il design
J  = 0.1;
T  = 810; %[N]
Cl=0.2*ones(size(el.r_bar));
for i=1:el.n_r
    [~,idx(i)]=min(abs(aero.Cl-Cl(i))); 
end
alpha=aero.alpha(idx);
el = el.Design(J,T,Cl,alpha,options);

%% MODEL 3D
data=importdata('NACA 16-212.dat');
x=data.data(:,1);
z=data.data(:,2);
figure
el.Model3D(x,z)



%% Function 
function CL=CL_(alpha)
global aero

    alpha=alpha*180/pi;
    if alpha <  aero.alpha(1)
        CL=aero.Cl(1);
    elseif alpha > aero.alpha(end)
        CL=aero.Cl(end);
    else 
        CL=interp1(aero.alpha,aero.Cl,alpha);
    end
end

function CD=CD_(alpha)
global aero
    alpha=alpha*180/pi;
    if alpha <  aero.alpha(1)
        CD=aero.Cd(1);
    elseif alpha > aero.alpha(end)
        CD=aero.Cd(end);
    else 
        CD=interp1(aero.alpha,aero.Cd,alpha);
    end
end