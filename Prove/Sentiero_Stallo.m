%% Sentiero di stallo
clc; clear; close all
rotore1 = Rotor();
% working conditions and other inputs
dim_vel = 50;
V_inf   = linspace(0.1,convvel(293,'km/h','m/s'),dim_vel);
Chi     = convang(5,'deg','rad');
f       = 3;
W       = 75278;
rotore1.theta_t = convang(-8,'deg','rad');
rotore1.h     = 0;
% properties
rotore1 = rotore1.r(linspace(0.1,1,100));
rotore1.R     = 7.6;
rotore1.N     = 3;
rotore1.c     = linspace(0.4,0.4,rotore1.n_r);
rotore1 = rotore1.ambient();
rotore1 = rotore1.mass_prop('G',8);
rotore1 = rotore1.rot_vel('omega',1,1);
% rotore1 = rotore1.BEMT_articulated('T',W,V_inf,Chi,f);
alpha_max_2D = convang(15,'deg','rad');

[s,r,c] = rotore1.sentiero_stallo(alpha_max_2D,'T',W,Chi,f);

%% POST PROCESSING
Video = VideoWriter('Sentiero_Stallo','Motion JPEG AVI');    % Open a Video file
Video_xy.Quality = 100;
open(Video);
writeVideo(Video,getframe(figure(1)))
V_inf = s.V_inf*[1.01 :0.02:1.2]';
s=alphamap(rotore1,'Solve',{'T',W,...
    V_inf,Chi,f,BEMTset_rotor()});
% plot di alpha_e_max
% Create polar data
[ra,psi] = meshgrid(rotore1.r_bar,s.options.Psi);
% Convert to Cartesian
x = ra.*cos(psi);   x1=0;
y = ra.*sin(psi);   y1=0;
for k= 1:length(V_inf)
    figure(k+1)
    hold on
    [r,c] = find(abs(s.alpha_e(:,:,k) - alpha_max_2D) < convang(0.1,'deg','rad'));
    for i=1:length(c)
        x1(i)=x(c(i),r(i));
        y1(i)=y(c(i),r(i));
    end
    plot(x1,y1,'--r')
    writeVideo(Video,getframe(figure(k+1)));
    pause(2)
    clear x1 y1
end
close(Video);


