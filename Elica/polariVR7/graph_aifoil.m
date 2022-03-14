% Rappresentazione profilo
clear all; close all; clc;

data = importdata('BOEING_VERTOL_VR-7.dat');
x    = data.data(:,1);
z    = data.data(:,2);

figure
plotta(x,z,{'$x/c$';'$z/c$'})
daspect([1 1 1])