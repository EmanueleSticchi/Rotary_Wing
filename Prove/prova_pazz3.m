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

syms w ap
syms sig lam1 lam2 sV_inf sOmR
assume(sV_inf>0)
assume(sOmR  >0)
assume(sig  >0)

eq1 = sig*lam1*((sV_inf + w)^2+(sOmR*(1-ap))^2) == ...
    4*(sV_inf + w)*w;
eq2 = sig*lam2*((sV_inf + w)^2+(sOmR*(1-ap))^2) == ...
    4*(sV_inf + w)*ap*sOmR;


eqns = [eq1, eq2];
vars = [w,ap];
[solv, solu] = solve(eqns,vars);

for jdx =1 :length(J)
    V_inf = J(jdx)*obj.n*obj.D;
    a  = 1;
    ap = 0;
    for idx = 1:length(obj.r_bar)
        res = 1;
        it = 0;
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
            
            OmR   = obj.omega*obj.r_bar(idx)*obj.R;
            wval  = subs(solv,[sig,lam1,lam2,sV_inf,sOmR],...
                [obj.sigma(idx),lambda1,lambda2,V_inf,OmR]);
            apval = subs(solu,[sig,lam1,lam2,sV_inf,sOmR],...
                [obj.sigma(idx),lambda1,lambda2,V_inf,OmR]);
            vW   = eval(wval);
            vap  = eval(apval);
            
            a_new  = vW(2)/V_inf;    err_a  = abs(a-a_new);
            ap_new = vap(1);         err_ap = abs(ap-ap_new);
            res    = max([err_ap,err_a]);
            a      = a_new;
            ap     = ap_new;
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