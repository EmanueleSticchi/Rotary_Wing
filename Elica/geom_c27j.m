%% Geometry Prop P92
clc; clear; close all
folder = 'immagini/Design/c27j/geom/';
load el_c27j.mat
theta_75 = convang(20,'deg','rad');
[~,idx] = min(abs(el.r_bar - 0.75));
el.theta = el.theta - el.theta(idx) + theta_75;

%% MODEL 3D
data=importdata('BOEING_VERTOL_VR-7.dat');
x=data.data(:,1);
z=data.data(:,2);
el.Model3D(x,z)
el.Model3D(x,z)
figure(1)
ax = gca;
view(ax,[-90 90]); 
f.Units = 'normalized';
f.Position = [0.13,0.11,0.775,0.815];
figure(2)
ax = gca;
view(ax,[25 30]); 
axis off

figure
plotta(el.r_bar,el.c,{'$ \bar{r}$';'c [m]'});
figure
plotta(el.r_bar,el.theta*180/pi,{'$ \bar{r}$';'$\theta$ [deg]'});
figure
plotta(el.r_bar,el.LAMBDA*180/pi,{'$ \bar{r}$';'$\Lambda$ [deg]'})
%% Save
count = 0;
for i =1:7
    count = count + 1;
    figure(i)
    FileName = sprintf(['geom','%d.eps'], count);
    ax = gca;
    exportgraphics(ax,[folder,FileName])
end

