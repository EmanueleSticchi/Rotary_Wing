%% Progetto turbina ( par 9.7 Tognaccini)
clc; clear; close all
load Aero_Du84-132V3_Re1e6.mat
%% Dati  iniziali -------------------------------------------------------
% geometria
obj.R     = 55;
obj.R_hub = 5;
obj.n_r   = 100;
obj.r     = linspace(obj.R_hub,obj.R,obj.n_r);
obj.r_bar = obj.r/obj.R;
obj.N     = 3;

% aerodinamica
obj.Cl    = @(alpha_rad) interp1(alpha,Cl,alpha_rad*pi/180); 

% funzionamento
V_inf      = 7;
lambda_opt = 6;
obj.omega  = lambda_opt*V_inf/obj.R;
OmR        = obj.omega*obj.r;
%% Design
Chi = OmR/V_inf;
obj.alpha_id = ones(obj.n_r,1)*convang(5,'deg','rad');
% Calettamento
toll = 1e-6;
phi0 = convang(50,'deg','rad');
phi  = zeros(obj.n_r,1);
for i=1:obj.n_r
     phi1 = phi0*0.7;
     it = 0;
     f0=func(phi0,Chi(i));
     if abs(f0)>toll
         f1=func(phi1,Chi(i));
         while abs(f1)> toll
             it = it+1;
             qk=(f1-f0)/(phi1-phi0);
             phi0=phi1; f0=f1;
             phi1=phi1-f1/qk;
             f1=func(phi1,Chi(i));
         end
     else
         phi1=phi0;
     end
     
     phi(i)       = phi1;
     % calettamento
     obj.theta(i) = phi(i) + obj.alpha_id(i);
     % corda
     obj.c(i) = (2*pi*V_inf)/(obj.N*obj.omega*obj.Cl(obj.alpha_id(i)))*(4*sin(phi(i))*(2*cos(phi(i))-1))/(1+2*cos(phi(i)));
end

%% Plotting 
figure
plot(obj.r,obj.c,'k')
xlabel('$r [m]$','Interpreter','latex','FontSize',12)
ylabel('$c [m]$','Interpreter','latex','FontSize',12,'Rotation',90);

figure
plot(obj.r,obj.theta*180/pi);
xlabel('$r [m]$','Interpreter','latex','FontSize',12)
ylabel('$\theta [deg]$','Interpreter','latex','FontSize',12,'Rotation',90);

save('turbina_proggettata.mat','obj')


%% Function
function f = func(phi,Chi)
% eq (9.33) pag 127
    f = Chi - (sin(phi)*(2*cos(phi)-1))/((1+2*cos(phi))*(1-cos(phi)));
end
