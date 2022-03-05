%%  Analisi di un'elica con la BEMT
clc; clear; close all

%%  DATI ----------------------------------------------------------
global aero
toll=1e-6;
obj=Elica();

% Geometrici
obj=obj.r_((0.1:0.01:1)');
obj.N=3;                         % Numero di pale    
obj.D=2;                         % Diametro elica, [m]
obj.R=0.5*obj.D;                  % Raggio dell'elica, [m]
obj.theta=pi/180*....            % Angolo di calettamento, [rad]
    linspace(50,15,obj.n_r)';
obj.c=0.1*ones(obj.n_r,1);        % Corda delle sezioni, [m]
obj=obj.sigma_();                 % Solidità
obj.LAMBDA = zeros(obj.n_r,1);
% Di funzionamento
obj=obj.rot_vel('RPM',2000);
obj= obj.altitude(0);
J=[0.5:0.05:1.3];
% Aerodinamici
% Cl=@(alpha) interp1(data1,data2,alpha);
% Cd=@(alpha) interp1(data1,data3,alpha);
load Aero_NACA16212.mat
obj.Cl=@(alpha,M,Re,r_bar) CL_(alpha);
obj.Cd=@(alpha,M,Re,r_bar) CD_(alpha);

%% Analisi BEMT --------------------------------------------------------
alpha0=-2*pi/180;
alpha1=10*pi/180;
options=BEMTset();
obj=obj.BEMT(J,alpha0,alpha1,options);
% options.Hub_correction='on';
% options.Cd_hub = 0.5;
% el=el.BEMT(J,alpha0,alpha1,options);

%% Post - Processing
s = obj.Analisi{1,1};

for jdx = 1:length(J)
   figure(1)
   plot(obj.r_bar,s.a)
   hold on
   figure(2)
   plot(obj.r_bar,s.ap)
   hold on
end





% an=2;
% for an=1:2
% figure(1)
% plotta(J,el.Analisi{an, 1}.CT,{'J = $ \frac{V_{\infty}}{nD}$';'$C_T$ = $\frac{T}{\rho n^2 D^4}$'})
% yline(0);
% 
% figure(2)
% plotta(J,el.Analisi{an, 1}.CP,{'J = $ \frac{V_{\infty}}{nD}$';'$C_P$ = $\frac{P}{\rho n^3 D^5}$'})
% yline(0);
% 
% figure(3)
% log=el.Analisi{an, 1}.CT>=0;
% plotta(J(log),el.Analisi{1, 1}.eta(log,1),{'J = $ \frac{V_{\infty}}{nD}$';'$\eta$ = $\frac{TV_{\infty}}{P}$'})
% yline(0);
% %
% idx=3;
% figure(4)
% plotta(el.r_bar,el.Analisi{an,1}.dCt_dr_bar(idx,:),{'$ \bar{r}$';'$\frac{dC_T}{d\bar{r}}$'})
% yline(0);
% figure(5)
% plotta(el.r_bar,el.Analisi{an,1}.dCq_dr_bar(idx,:),{'$ \bar{r}$';'$\frac{dC_Q}{d\bar{r}}$'})
% yline(0);
% figure(6)
% plotta(el.r_bar,el.Analisi{an,1}.dCp_dr_bar(idx,:),{'$ \bar{r}$';'$\frac{dC_P}{d\bar{r}}$'})
% yline(0);
% figure(7)
% plotta(el.r_bar,el.Analisi{an,1}.eta_e(idx,:),{'$ \bar{r}$';'$\eta_e$'})
% yline(0);
% end


%% Function ---------------------------------------------------------
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