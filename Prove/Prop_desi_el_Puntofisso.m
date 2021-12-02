%%                      Propellers Design, T assigned
clc; clear; close all
global aero
%%  DATA ------------------------------------------------------------------
el=Elica();

% Geometry
el=el.r_((0.1:0.01:1)');
el.N=4;                         % Blade Number    
el.D=1.66;                         % Propeller diameter, [m]
el.R=0.5*el.D;                  % Propeller radius, [m]

% Aerodinamici
% Cl=@(alpha) interp1(data1,data2,alpha);
% Cd=@(alpha) interp1(data1,data3,alpha);
load Aero_NACA16212.mat
el.Cl=@(alpha) CL_(alpha);
el.Cd=@(alpha) CD_(alpha);

% Operating condition
el=el.rot_vel('RPM',2000);
el.rho=1.23;
J=0.;
V_inf=J*el.n*el.D;
T=1000;                        % Thrust, [N]
Cl=0.2*ones(size(el.r_bar));
for i=1:el.n_r; [~,idx(i)]=min(abs(aero.Cl-Cl(i))); end
alpha=aero.alpha(idx);
% alpha=0.5*Cl/pi*pi/180;

%% Design
el=el.Design(J,T,Cl,alpha);

%% Analisys
alpha0=-20*pi/180;
alpha1=10*pi/180;
J=[0.:0.1:1.5];
options=BEMTset();
options.P_correction='on';
el=el.BEMT(J,alpha0,alpha1,options);
%% PLotting 
data=importdata('NACA 16-212.dat');
x=data.data(:,1);
z=data.data(:,2);
figure
plot(x,z)
daspect([1 1 1])
figure
el.Model3D(x,z)
%
an=1;
figure
plotta(J,el.Analisi{an, 1}.CT,{'J = $ \frac{V_{\infty}}{nD}$';'$C_T$ = $\frac{T}{\rho n^2 D^4}$'})
yline(0);

figure
plotta(J,el.Analisi{an, 1}.CP,{'J = $ \frac{V_{\infty}}{nD}$';'$C_P$ = $\frac{P}{\rho n^3 D^5}$'})
yline(0);

figure
log=el.Analisi{an, 1}.CT>=0;
plotta(J(log),el.Analisi{1, 1}.eta(log,1),{'J = $ \frac{V_{\infty}}{nD}$';'$\eta$ = $\frac{TV_{\infty}}{P}$'})
yline(0);




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

function M=R2(alpha)
M=[cos(alpha)  0  -sin(alpha);...
    0          1        0
   sin(alpha)  0    cos(alpha) ];
end
function M=R3(alpha)
M=[ cos(alpha) sin(alpha) 0;...
    -sin(alpha) cos(alpha) 0;...
        0          0     1];
end