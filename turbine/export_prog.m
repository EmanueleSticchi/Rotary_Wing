%% Export turbina for Qblade in Aerodyn format
clc; clear; close all
%
load turbina_proggettata.mat
folder   = 'turb_prog_aerodyn/'; 
nomefile = 'aerodyn.ipt';
id       = fopen([folder,nomefile], 'w+');
dr       = obj.r(2:end) - obj.r(1:end-1);
dr       = [ dr(1), dr];
mat = [obj.r;obj.theta*180/pi;dr;obj.c];
fprintf(id, '%f  %f   %f   %f   1   NOPRINT\n',mat);
fclose(id)