%%                      Propellers Design, T assigned
clc; clear; close all
global aero
%%  DATA ------------------------------------------------------------------
toll=1e-6;
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
J=0.5;
V_inf=J*el.n*el.D;
T=1000;                        % Thrust, [N]
Cl=0.2*ones(size(el.r_bar));
for i=1:el.n_r; [~,idx(i)]=min(abs(aero.Cl-Cl(i))); end
alpha=aero.alpha(idx);
% alpha=0.5*Cl/pi*pi/180;

%% Design
el=el.Design(J,T,Cl,alpha);

%% Analisys
alpha0=-2*pi/180;
alpha1=10*pi/180;
J=[0.2:0.05:2.5];
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



% for i=1:el.n_r
%     M=R2(el.theta(i));
%     M=M([1 3],[1 3]);
%     % set pitch
%     for j=1:length(x)
%         data_rot(j,:)=M*data.data(j,:)';
%     end
%     % scale airfoil
%     x_rot=data_rot(:,1)*el.c(i);
%     z_rot=data_rot(:,2)*el.c(i);
%     X(:,i)=x_rot-mean(x_rot);
%     Z(:,i)=z_rot;
% end
% 
% Y=el.r_bar'.*el.R;
% Y=repmat(Y,length(x),1);
% figure
% light('Style','local','Position',[1 -1 0]);
% 
% s=surf(X,Y,Z,'FaceColor',[0.65 0.65 0.65],'FaceLighting','gouraud','EdgeColor','none');
% daspect([1 1 1])
% camlight right 
% % material dull
% 
% % Hub(cilinder)
% m=100;
% % create Hub disc
% xh=linspace(-el.r_bar(1),el.r_bar(1),m)*el.R;
% yh=sqrt((el.r_bar(1)*el.R)^2-xh.^2);
% xh=[xh,flip(xh)]';
% yh=[yh,-yh]';
% Xh=repmat(xh,1,3);
% Yh=repmat(yh,1,3);
% % estrusion of the th Hub Disc
% Zh=repmat(-0.1:0.1:0.1,length(xh),1);
% hold on
% surf(Xh,Yh,Zh)
% daspect([1 1 1])
% camlight right 
% % Hub(parabolic)
% m=100;
% % create Hub disc
% xh=linspace(-el.r_bar(1),el.r_bar(1),m)'*el.R;
% yh=xh;
% % xh=[xh,flip(xh)]';
% % yh=[yh,-yh]';
% % Xh=repmat(xh,1,3);
% % Yh=repmat(yh,1,3);
% [Xh,Yh]=meshgrid(xh,yh);
% % estrusion of the th Hub Disc
% k=20;
% Zh=-k*Xh.^2-k*Yh.^2+0.2;
% hold on
% surf(Xh,Yh,Zh)
% daspect([1 1 1])
% 
% 
% 
% 
% % other blade
% ang_blade=2*pi/el.N;
% ang=0;
% 
% for i =1:el.N
%     ang=ang+ang_blade;
%     xyz=[X(:),Y(:),Z(:)];
%     xyz_rot=xyz*R3(ang);
%     X_rot=reshape(xyz_rot(:,1),[length(x),el.n_r]);
%     Y_rot=reshape(xyz_rot(:,2),[length(x),el.n_r]);
%     Z_rot=reshape(xyz_rot(:,3),[length(x),el.n_r]);
%     s=surf(X_rot,Y_rot,Z_rot,'FaceColor',[0.65 0.65 0.65],'FaceLighting','gouraud','EdgeColor','none');
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