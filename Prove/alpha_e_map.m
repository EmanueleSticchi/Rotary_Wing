clear all; close all; clc;

rotore1 = Rotor();
% working conditions and other inputs
dim_vel = 50;
V_inf   = linspace(0.1,convvel(293,'km/h','m/s'),dim_vel);
Chi     = convang(5,'deg','rad');
f       = 3;
W       = 75278;
theta_t = convang(-8,'deg','rad');
rotore1.h     = 0;
% properties
rotore1 = rotore1.r(linspace(0.1,1,100));
rotore1.R     = 7.6;
rotore1.N     = 3;
rotore1.c     = linspace(0.4,0.4,rotore1.n_r);
rotore1.theta = pi/180*linspace(13.3,9,rotore1.n_r);
% obj.I = convmass(obj.Ixx,'slug','kg');
% obj.I = convlength(convlength(obj.Ixx,'ft','m'),'ft','m');
% function recall
rotore1 = rotore1.ambient();

rotore1 = rotore1.mass_prop('G',8);
rotore1 = rotore1.rot_vel('omega',1,1);
rotore1 = rotore1.BEMT_articulated(V_inf,Chi,f,W,theta_t);
% W1 = W*1.3;
% rotore1 = rotore1.BEMT_articulated(V_inf,Chi,f,W1,theta_t);
% W2 = W*1.5;
% rotore1 = rotore1.BEMT_articulated(V_inf,Chi,f,W2,theta_t);

%% Controllo
% s = rotore1.Analisi_articulated{1,1};
% Psi=s.options.Psi;
% alpha_e=zeros(rotore1.n_r,length(Psi),length(s.lam_Vec));
% idxV = 16;
% b     =  s.beta0_Vec(idxV) + ...
%     s.beta1c_Vec(idxV)*cos(Psi) +...
%     s.beta1s_Vec(idxV)*sin(Psi);
% 
% b_dot = -s.beta1c_Vec(idxV)*sin(Psi) +...
%     s.beta1s_Vec(idxV)*cos(Psi);
% figure
% plot(Psi*180/pi,b*180/pi,'b',Psi*180/pi,b_dot*180/pi,'r')
% figure
% plot(rotore1.r_bar,s.theta(:,idxV)*180/pi)
% figure
% plot(s.mu,s.lam_Vec,s.mu(idxV),s.lam_Vec(idxV),'*','MarkerSize',10)
% figure
% i=r;
% j=c;
% alpha_e(i,j,idxV) = s.theta(i,idxV) -...
%             (s.lam_Vec(idxV)  +...
%             b_dot(j)*rotore1.r_bar(i)/rotore1.omega+...
%             b(j)*s.mu(idxV)*cos(Psi(j)))/(...
%             rotore1.r_bar(i) + s.mu(idxV)*sin(Psi(j)))




% for i=1:rotore1.n_r
%     for j=1:length(Psi)
%         alpha_e(i,j,idxV) = s.theta(i,idxV) -...
%             (s.lam_Vec(idxV)  +...
%             b_dot(j)*rotore1.r_bar(i)/rotore1.omega+...
%             b(j)*s.mu(idxV)*cos(Psi(j)))/(...
%             rotore1.r_bar(i) + s.mu(idxV)*sin(Psi(j)));
%     end
% end





%% Alpha_e map 
s = rotore1.Analisi_articulated{1,1};
idxV = 50;
figure
alpha_e=s.alpha_e(:,:,idxV)*180/pi;
% alpha_e(alpha_e < -10) = -10;
% contour(rotore1.r_bar,s.options.Psi,alpha_e')

% Create polar data
[r,psi] = meshgrid(rotore1.r_bar,s.options.Psi);
% Convert to Cartesian
x = r.*cos(psi);
y = r.*sin(psi);

% define polar axes
h = polar(x,y);
hold on;
polar(s.options.Psi,rotore1.r_bar(1)*ones(length(s.options.Psi),1),'k')
polar(s.options.Psi,rotore1.r_bar(end)*ones(length(s.options.Psi),1),'k')
% contourf(x,y,alpha_e');
pc= pcolor(x,y,alpha_e');
contour(x,y,alpha_e','k','ShowText','on');

shading interp
% colormap 'hsv'
cbar=colorbar(gca);
cbar.Label.String = '\alpha_e';
cbar.Label.FontSize= 16;
% cbar.Limits = [-10 10];



% Hide the POLAR function data and leave annotations
set(h,'Visible','off')
% Turn off axes and set square aspect ratio
axis off
axis image
view([90 90])
title(['\alpha_{e_{max}} = ',num2str(max(alpha_e,[],'all')),' deg'])

m=max(abs(alpha_e),[],'all');
[r,c] = find(abs(alpha_e)==m)


%% Graphics
% figure;
% plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{1,1}.Pc_Vec,'-k');
% hold on
% plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{2,1}.Pc_Vec,':k');
% plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{3,1}.Pc_Vec,'.-k');
% 
% figure;
% plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{1,1}.beta0_Vec,'-k');
% hold on
% plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{1,1}.beta1c_Vec,':k');
% plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{1,1}.beta1s_Vec,'.-k');
% plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{2,1}.beta0_Vec,'-r');
% plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{2,1}.beta1c_Vec,':r');
% plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{2,1}.beta1s_Vec,'.-r');
% plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{3,1}.beta0_Vec,'-g');
% plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{3,1}.beta1c_Vec,':g');
% plot(V_inf/(rotore1.R*rotore1.omega),rotore1.Analisi_articulated{3,1}.beta1s_Vec,'.-g');

