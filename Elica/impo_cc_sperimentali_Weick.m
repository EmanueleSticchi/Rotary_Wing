%% importa curve caratteristiche sperimentali
clc;clear; close all
data = importdata('CC_sperimentali.txt');

CC_exp.J   = data(1:15);
CC_exp.Ct  = data(16:30);
CC_exp.Cp  = data(31:45);
CC_exp.eta = data(46:end);

save('CC_experimential_Weick.mat','CC_exp');