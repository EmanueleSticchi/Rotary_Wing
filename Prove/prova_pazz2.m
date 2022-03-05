%%  Analisi di un'elica con la BEMT
clc; clear; close all

%%  DATI ----------------------------------------------------------
global aero
toll=1e-6;
obj=Elica2();

% Geometrici
obj=obj.r_((0.1:0.01:1)');
obj.N=3;                         % Numero di pale    
obj.D=2;                         % Diametro elica, [m]
obj.R=0.5*obj.D;                  % Raggio dell'elica, [m]
obj.theta=pi/180*....            % Angolo di calettamento, [rad]
    linspace(50,15,obj.n_r)';
obj.c=0.1*ones(obj.n_r,1);        % Corda delle sezioni, [m]
obj=obj.sigma_();                 % Solidità

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
%% Analisi con BEMT ( nuova procedura) 
toll = 1e-6;
for jdx = 1: length(J)
    V_inf = J(jdx)*obj.n*obj.D;
    % valori di primo tentativo per a e ap
    a  = 1;
    ap = 1;
    for idx = 1:obj.n_r
        res = 1;
        it = 0;
%         a  = 0.01*V_inf;
%         ap = 0;
        while res > toll
            it = it+1;
            
            phi     = atan2(V_inf*(1+a),(obj.omega*obj.r_bar(idx)*obj.R*(1 - ap)));
            V_eff   = sqrt((V_inf*(1+a))^2 + (obj.omega*obj.r_bar(idx)*obj.R*(1 - ap))^2);
            M       = V_eff/obj.sound_vel;
            Re      = V_eff * obj.rho * obj.c(idx)/ obj.mu_visc;
            alpha   = obj.theta(idx) - phi;
            lambda1 = (obj.Cl(alpha,M,Re,obj.r_bar(idx))*cos(phi)...
                -obj.Cd(alpha,M,Re,obj.r_bar(idx))*sin(phi));
            lambda2 = (obj.Cl(alpha,M,Re,obj.r_bar(idx))*sin(phi)...
                +obj.Cd(alpha,M,Re,obj.r_bar(idx))*cos(phi));
% %             metodo con a e ap classico (non converge)(non valida a punto
% %             fisso)
%             ka      = 0.25*obj.sigma(idx)*lambda1/sin(phi)^2;
%             kap     = 0.5*obj.sigma(idx)*lambda2/sin(2*phi);
%             
%             ap_new  = kap/(1+kap);  err_ap = abs(ap-ap_new);
%             a_new   = ka/(1-ka);    err_a  = abs(a-a_new);
%             res     = max([err_a,err_ap]);

            % metodo con a e ap soluzioni del sistema (valida a punto fisso)
            OmR = obj.omega*obj.r_bar(idx)*obj.R;
            A   = obj.sigma(idx)*(lambda1 + lambda2^2/lambda1) - 4;
            B   = 2*(obj.sigma(idx)*(lambda1*V_inf - lambda2*OmR) - 2*V_inf);
            C   = obj.sigma(idx)*(lambda1*V_inf^2 + lambda1*OmR^2);
            wp   = (-B + sqrt(B^2 -4*A*C))/(2*A);
            wm   = (-B - sqrt(B^2 -4*A*C))/(2*A);
            a_new  = wp/V_inf;                err_a  = abs(a-a_new);
            ap_new = wm/OmR*lambda2/lambda1;  err_ap = abs(ap-ap_new);
            res     = max([err_a,err_ap]);
%             res = err_a;

            a       = a_new;
            ap      = ap_new;
        end
        clc
        disp(['res = ',num2str(res)])
        disp(['it = ',num2str(it)])
        disp(['a = ',num2str(a)])
        disp(['ap = ',num2str(ap)])
        % save for plot
        a_vec(idx)  = a;
        ap_vec(idx) = ap;
        it_vec(idx) = it;
        re_vec(idx) = res;
    end
    figure(1)
    plot(obj.r_bar,a_vec)
    hold on
    figure(2)
    plot(obj.r_bar,ap_vec)
    hold on
    figure(3)
    plot(obj.r_bar,it_vec)
    hold on
    figure(4)
    plot(obj.r_bar,re_vec)
    hold on
end






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