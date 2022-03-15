classdef Elica
    properties  
    % ---------------------------------------------------------------------
    % Geometria
    % ---------------------------------------------------------------------
        N     {mustBeInteger, mustBeFinite}        % Numero di pale
        D     {mustBePositive, mustBeFinite}       % Diametro elica, [m]
        A_D   {mustBePositive, mustBeFinite}       % Area disco elica,[m2]
        R     {mustBePositive, mustBeFinite}       % Raggio dell'elica, [m]                
        theta (:,1){mustBeReal, mustBeFinite}%Angoli di calettamento, [rad]
        c     (:,1){mustBeNonnegative, mustBeFinite}% Corda delle sezioni, [m]
        sigma (:,1){mustBeNonnegative, mustBeFinite}       % Solidità
        LAMBDA (:,1){mustBeReal, mustBeFinite}      % Angolo di Freccia [rad]
    % --------------------------------------------------------------------- 
    % Funzionamento
    % ---------------------------------------------------------------------
        RPM   {mustBePositive, mustBeFinite}       % Giri al minuto
        n     {mustBePositive, mustBeFinite}       % Giri al secondo
        omega {mustBePositive, mustBeFinite}%Velocità di rotazione, [rad/s]
        h {mustBeNonnegative, mustBeFinite} = 0 % Quota di funzionamento [m]
    % ---------------------------------------------------------------------
    % Aerodinamica
    % ---------------------------------------------------------------------
        Cl= @(alpha,r_bar,M,Re) 2*pi*alpha;
        Cd= @(alpha,r_bar,M,Re) 0.01*alpha./alpha;
        
    end
    properties(SetAccess = private,GetAccess=public)
        % Vettore delle stazioni radiali
        r_bar (:,1) {mustBeInRange(r_bar,0,1,'exclude-lower')} %exclusive
        n_r   
        %
        n_analisi = 0
        % ambient conditons
        rho       {mustBePositive, mustBeFinite}
        press     {mustBePositive, mustBeFinite}
        sound_vel {mustBePositive, mustBeFinite}
        temp      {mustBePositive, mustBeFinite}
        mu_visc   {mustBePositive, mustBeFinite}
        % funzionamento
        OmR
        % ---------------------------------------------------------------------
        % Analisi
        % ---------------------------------------------------------------------
        Analisi
    end
    properties(SetAccess=private,GetAccess=private)
        Des
    end
% ---------------------------------------------------------------------
    methods
        % imposta il dominio radiale
        function obj = r_(obj,vec_r)
            obj.r_bar=vec_r;
            obj.n_r=length(obj.r_bar);
        end
        % Compute some mass and geometric property. This function needs to
        % be called before doing any other calculations but still after the
        % definition of the main properties(N,c,r_bar,R,omega). The function needs to be
        % called only once. The function computes any derived property.
        function obj = derived_properties(obj)
            obj.D     = obj.R*2;
            obj.A_D   = pi*obj.R*obj.R;
            if ~isempty(obj.c)
                obj.sigma = obj.N/(2*pi)*obj.c.*obj.r_bar.^-1*obj.R;
            end
        end
        % calcolo delle velocità di rotazione
        function obj = rot_vel(obj,valIN,val)
            switch valIN
                case 'RPM'
                    obj.RPM   = val;
                    obj.n     = obj.RPM/60;
                    obj.omega = obj.n*2*pi;
                case 'n'
                    obj.n     = val;
                    obj.RPM   = obj.n*60;
                    obj.omega = obj.n*2*pi;
                case 'omega'
                    obj.omega = val;
                    obj.n     = obj.omega/(2*pi);
                    obj.RPM   = obj.n*60;
                otherwise 
                    mustBeMember(valIN,{'RPM','n','omega'})
            end
            obj.OmR   = obj.omega*obj.r_bar*obj.R;
        end
        % assegna le prestazioni aerodinamiche di profilo in funzione di 
        % dati tabulari... Da modificare
        function obj = set_aero(obj,valpha,vCl,vCd)
            % v_alpha   vettore degli angoli di attacco in rad
            % vCl       vettore dei Cl (stessa dimensione di v_alpha)
            % vCd       vettore dei Cd (stessa dimensione di v_alpha)
            obj.Cl=@(alpha) interp1(valpha,vCl,alpha);
            obj.Cd=@(alpha) interp1(valpha,vCd,alpha);
        end
        % -----------------------------------------------------------------
        % Ambient conditions
        % -----------------------------------------------------------------
        function obj = altitude(obj,h)
            obj.h=h;
            [obj.temp, obj.sound_vel, obj.press, obj.rho] = atmosisa(obj.h);
            T0  = 288.15;
            mu0 = 1.79e-5;
            obj.mu_visc=mu0*(obj.temp/T0)^1.5*((T0+110)/(obj.temp+110));
        end
        
        % ----------------------------------------------------------------
        %% Teoria dell'elemento di pala generale
        % ----------------------------------------------------------------
        function [f,J_,alpha,phi,a,ap,lambda1,lambda2,M,Re]=func(obj,alpha,J,idx,options)
        % Funzione da annullare per il calcolo delle prestazioni dell'elica
        % Input:
        % - alpha,           valore di alpha di tentativo;
        % - J,               Rapporto di funzionamento in corrispondenza
        %                    del quale si vogliono calcolare le prestazioni
        % - idx,             Indice della stazione radiale
        % - options          Analisys options (see BEMTset)
        % -----------------------------------------------------------------
            V_inf = J *obj.n*obj.D;
            V_eff = sqrt(V_inf^2 + (obj.omega*obj.r_bar(idx)*obj.R)^2)*cos(obj.LAMBDA(idx)); % trascuro a e a' rispetto ad 1
            M  = V_eff/obj.sound_vel;
            Re = V_eff * obj.rho * obj.c(idx)/ obj.mu_visc;
            % calcolo dell'angolo di Inflow
            phi=obj.theta(idx)-alpha;
            % calcolo dei coefficienti aerodinamici
            Cl = obj.Cl(alpha,obj.r_bar(idx),M,Re);
            Cd = obj.Cd(alpha,obj.r_bar(idx),M,Re);
            
            lambda1 = Cl*cos(phi) - Cd*sin(phi);
            lambda2 = Cl*sin(phi) + Cd*cos(phi);
            % calcolo delle induzioni
%             ka  = 0.25*obj.sigma(idx)*lambda1/(sin(phi)^2);
            ka  = 0.25*obj.sigma(idx)*lambda1/(1- cos(phi)^2);
            kap = 0.5*obj.sigma(idx)*lambda2/sin(2*phi);

            ap = kap/(1+kap);
            a  = ka/(1-ka);
            
            % calcolo del rapporto di avanzamento
            J_=pi*obj.r_bar(idx)*(1-ap)*tan(phi)/(1+a);
            % funzione da annullare
            f=J-J_;
        end
        function [f,J_,alpha,phi,a,ap,lambda1,lambda2,M,Re]=BEMT_rJ_fix(obj,...
                                                alpha0,alpha1,J,idx,options)
        % Calcola le prestazioni dell'elica per una sola stazione radiale
        %  e per un solo valore del rapporto di funzionamento con il metodo
        %  delle Secanti.
        % Input:
        % - alpha0 e alpha1       Valori di tentativo iniziali dell'angolo 
        %                           d'attacco; 
        % - J,                    Rapporto di funzionamento in 
        %                           corrispondenza del quale si vogliono  
        %                           calcolare le prestazioni;
        % - idx,                  Indice della stazione radiale
        % - options,              Analisys options (see BEMTset)
        % -----------------------------------------------------------------    
            it = 0;
            f0=obj.func(alpha0,J,idx,options);
            if abs(f0)>options.toll
                f1=obj.func(alpha1,J,idx,options);
                while abs(f1)>options.toll
                    it = it+1;
                    qk=(f1-f0)/(alpha1-alpha0);
                    alpha0=alpha1; f0=f1;
                    alpha1=alpha1-f1/qk;
                    f1=obj.func(alpha1,J,idx,options);
                end
            else
                alpha1=alpha0;
            end
            [f,J_,alpha,phi,a,ap,lambda1,lambda2,M,Re]=obj.func(alpha1,J,idx,options);

        end
        function obj=BEMT(obj,J,alpha0,options)
        % Calcola le prestazioni dell'elica per tutte le stazioni radiali
        %  e per una serie di valori del rapporto di funzionamento, con il 
        %  metodo delle Secanti.
        % Input:
        % - alpha0 e alpha1       Valori di tentativo iniziali dell'angolo 
        %                           d'attacco; [rad] (default: 0 e 10°)
        % - J,                    Valori del Rapporto di funzionamento in 
        %                           corrispondenza dei quali si vogliono  
        %                           calcolare le prestazioni;
        % - options               Analisys options (see BEMTset)
        % -----------------------------------------------------------------
        arguments
            obj
            J      (:,1){mustBeNonnegative,mustBeFinite}
            alpha0 (1,1){mustBeReal,mustBeFinite}=0
            options = BEMTset();
        end
        % -----------------------------------------------------------------
        alpha  = alpha0;
            for jdx=1:length(J)
                if jdx>1
                    alpha = s.alpha(jdx-1,1);
                end
                for idx =1:obj.n_r
                    alpha1 = 1.3*alpha;
                    [~,~,alpha,phi,a,ap,lambda1,lambda2,M,Re]=BEMT_rJ_fix(obj,...
                                           alpha,alpha1,J(jdx),idx,options);
                    s.alpha(jdx,idx)   = alpha;
                    s.phi(jdx,idx)     = phi;
                    s.a(jdx,idx)       = a;
                    s.ap(jdx,idx)      = ap;
                    s.lambda1(jdx,idx) = lambda1;
                    s.lambda2(jdx,idx) = lambda2;
                    s.Mach(jdx,idx)    = M;
                    s.Re(jdx,idx)      = Re;
                    
                    % calcolo delle prestazioni
                    s.dCt_dr_bar(jdx,idx) = pi/4*lambda1*obj.sigma(idx)*obj.r_bar(idx)*...
                        ((J(jdx)*(1+a))^2 + (pi*obj.r_bar(idx)*(1-ap))^2);
                    s.dCq_dr_bar(jdx,idx) = pi/8*lambda2*obj.sigma(idx)*obj.r_bar(idx)^2*...
                        ((J(jdx)*(1+a))^2 + (pi*obj.r_bar(idx)*(1-ap))^2);
%                     s.dCt_dr_bar(jdx,idx) = pi^3/4*obj.sigma(idx)*...
%                         lambda1*obj.r_bar(idx)^3*(1-ap)^2/(cos(phi)^2);
%                     s.dCq_dr_bar(jdx,idx) = pi^3/8*obj.sigma(idx)*...
%                         lambda2*obj.r_bar(idx)^4*(1-ap)^2/(cos(phi)^2);                   
                    s.dCp_dr_bar(jdx,idx) = 2*pi*s.dCq_dr_bar(jdx,idx);
                    s.eta_e(jdx,idx)      = J(jdx)/(pi*obj.r_bar(idx))*lambda1/lambda2;
%                     s.eta_e(jdx,idx)      = (1-ap)/(1+a)*lambda1/lambda2*tan(phi); 
                end 
                if isequal(options.P_correction,'on')
                    F                   = obj.F_(J(jdx)/pi);
                    s.dCt_dr_bar(jdx,:) = s.dCt_dr_bar(jdx,:)'.*F;
                    s.dCq_dr_bar(jdx,:) = s.dCq_dr_bar(jdx,:)'.*F;
                    s.dCp_dr_bar(jdx,:) = s.dCp_dr_bar(jdx,:)'.*F;
                end
                s.CT(jdx,1)=obj.simpsons(s.dCt_dr_bar(jdx,:),obj.r_bar(1),obj.r_bar(end));
                if isequal(options.Hub_correction,'on')
                    s.DCT(jdx,1) = -pi/8*(obj.r_bar(1))^2*...
                                J(jdx)^2*options.Cd_hub;
                    s.CT(jdx,1)  = s.CT(jdx,1) + s.DCT(jdx,1);
                end
                s.CQ(jdx,1)=obj.simpsons(s.dCq_dr_bar(jdx,:),obj.r_bar(1),obj.r_bar(end));
                s.CP(jdx,1)=2*pi*s.CQ(jdx,1);
                s.eta(jdx,1)=s.CT(jdx,1)/s.CP(jdx,1)*J(jdx);
            end
            s.J=J;
            s.options=options;
            if isequal(options.Design,'on')
                obj.Des=s.CT;
            else
                obj.n_analisi=obj.n_analisi+1;
                obj.Analisi{obj.n_analisi,1}=s;
            end
        end
        function I = simpsons(obj,f,a,b)
            h=(b-a)/(length(f)-1);
            I= h/3*(f(1)+2*sum(f(3:2:end-2))+4*sum(f(2:2:end-1))+f(end));
        end
        % -----------------------------------------------------------------
        %% Propeller Design
        % -----------------------------------------------------------------
        function obj=Design(obj,J,T,Cl,alpha,options)
        % Design of a propeller that minimize the required power in certain 
        %  operating condition (fixed J and Thrust T).
        % INPUT:
        % - J                   Design advance ratio
        % - T                   Design Trust [N] 
        % - Cl                  Design lift coefficient for each radial 
        %                           station
        % - alpha               Angle of attack corrisponding to Cl
        % - options             Analisys options (see BEMTset)
        % -----------------------------------------------------------------
            arguments
                obj
                J  {mustBeNonnegative,mustBeFinite}
                T  {mustBePositive,mustBeFinite}
                Cl {mustBeReal,mustBeFinite}
                alpha {mustBeReal,mustBeFinite}
                options=BEMTset();
            end
            options.Design='on';
            options.P_correction='on';
            if isempty(obj.LAMBDA)
                obj.LAMBDA = zeros(obj.n_r,1); 
            end
        % -----------------------------------------------------------------
            if isequal(options.Freccia_opt,'on')
                for i =1:obj.n_r
                    M_eff = sqrt((obj.OmR(i))^2+(J*obj.n*obj.D)^2)/obj.sound_vel;
                    if M_eff > options.M_lim
                        obj.LAMBDA(i) = acos(options.M_lim/M_eff);
                    else
                        obj.LAMBDA(i) = 0;
                    end
                end
            end
            V_inf=J*obj.n*obj.D;
            A=pi*obj.R^2;
            CT= T/(obj.rho*obj.n^2*obj.D^4);
            % X=el.omega*el.R/V_inf*el.r_bar;
            if V_inf < 1e-4
                w0=sqrt(T/(2*obj.rho*A));
            else
                a=-0.5*(1-sqrt(1+2*T/(obj.rho*V_inf^2*A)));
                w0=a*V_inf;        
            end
            f0=obj.funcDes(w0,V_inf,Cl,alpha,CT,options);
            if abs(f0)>options.toll
                w01=1.5*w0;
                f1=obj.funcDes(w01,V_inf,Cl,alpha,CT,options);
                while abs(f1)>options.toll
                    qk=(f1-f0)/(w01-w0);
                    w0=w01; f0=f1;
                    w01=w01-f1/qk;
                    f1=obj.funcDes(w01,V_inf,Cl,alpha,CT,options);
                end
            else
                w01=w0;
            end
            [~,obj]=obj.funcDes(w01,V_inf,Cl,alpha,CT,options);
        end
        function [f,obj]=funcDes(obj,w0,V_inf,Cl_id,alpha_id,CT,options)
            phi = atan2((V_inf+w0),(obj.omega*obj.r_bar*obj.R));
            if V_inf< 1e-4
                OR_wo=obj.omega*obj.r_bar*obj.R/w0;
                w=w0*(1+OR_wo.^-2).^-1;
                ap=1*(1+OR_wo.^2).^-1;
                Ve  =sqrt((obj.omega*obj.r_bar*obj.R).^2.*(1-ap).^2+w.^2);
            else
                Chi=obj.omega*obj.r_bar*obj.R/V_inf;
                a   = w0/V_inf*Chi.^2.*((1+w0/V_inf)^2+Chi.^2).^-1;
                ap  = w0/V_inf*(1+w0/V_inf)^2*((1+w0/V_inf)^2+Chi.^2).^-1;
                Ve  =sqrt((obj.omega*obj.r_bar*obj.R).^2.*(1-ap).^2+V_inf^2*(1+a).^2);
            end
            Gamma = obj.omega*(obj.r_bar*obj.R).^2.*(4*pi*obj.r_bar.^2.*ap)...
                .*obj.F_(V_inf/obj.omega/obj.R);
            obj.sigma = 1/(pi*obj.R)*(Ve.*obj.r_bar.*Cl_id).^-1.*Gamma;
            obj.c=(2*pi)*obj.sigma.*obj.r_bar*obj.R/obj.N;
            obj.theta=phi+alpha_id;

            obj=obj.BEMT(V_inf/(obj.n*obj.D),-2*pi/180,options);
            f=obj.Des/CT -1;
        end
        function F=F_(obj,lambda)
            % Prandtl Correction F (computed for all r_bar)
            % Input:
            % - lambda                     advance ratio (see pag. 44)[1,1]
            % Output :
            % - F                          Prandtl correction factor
            %                               [obj.n_r,1]
            F=2/pi*acos(exp(0.5*obj.N/lambda*(obj.r_bar-1)));
            if lambda ==0
                F(isnan(F))=1;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PLOTTING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Model3D(obj,x,z) %,center
            % Plot 3D model of the propeller
            % INPUT:
            % - x e z               Airfoil coordinate (in percent of chord)
            % - center              x coordinate of the center of the
            %                       propeller (by default is assumed the 
            %                       center of the local chord)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             arguments
%                 obj
%                 x
%                 z
%                 center
%             end
            %--------------------------------------------------------------
            r_star = obj.r_bar(obj.LAMBDA == 0);
            r_star = r_star(end);       % individuo la stazione dove parte la freccia
            % plot the first blade
            xz=[x,z];
            for i=1:obj.n_r
                
                M=R2(obj.theta(i));
                M=M([1 3],[1 3]);
                % set pitch
                for j=1:length(x)
                    data_rot(j,:)=M*xz(j,:)';
                end
                % scale airfoil
                x_rot  = data_rot(:,1)*obj.c(i);
                z_rot  = data_rot(:,2)*obj.c(i);
                X(:,i) = x_rot - ... 
                    mean(x_rot) -...
                        (obj.r_bar(i) - r_star)*obj.R*tan(obj.LAMBDA(i))*cos(obj.theta(i));
                Z(:,i)=z_rot - (obj.r_bar(i) - r_star)*obj.R*tan(obj.LAMBDA(i))*sin(obj.theta(i));
            end

            Y=obj.r_bar'.*obj.R;
            Y=repmat(Y,length(x),1);
            % plot di una sola pala
            figure
            light('Style','local','Position',[1 -1 0]);
            s=surf(X,Y,Z,'FaceColor',[0.65 0.65 0.65],'FaceLighting',...
                'gouraud','EdgeColor','none');
            daspect([1 1 1])
            camlight right 
            % material dull
            xlabel('X [m]')
            ylabel('Y [m]')
            % plot di tutte le pale
            figure
            light('Style','local','Position',[1 -1 0]);
            ang_blade=2*pi/obj.N;
            xyz=[X(:),Y(:),Z(:)];
            ang = 0;
            hold on
            grid on
            for i =1:obj.N
                xyz_rot=xyz*R3(ang);
                X_rot=reshape(xyz_rot(:,1),[length(x),obj.n_r]);
                Y_rot=reshape(xyz_rot(:,2),[length(x),obj.n_r]);
                Z_rot=reshape(xyz_rot(:,3),[length(x),obj.n_r]);
                s=surf(X_rot,Y_rot,Z_rot,'FaceColor',[0.65 0.65 0.65],...
                    'FaceLighting','gouraud','EdgeColor','none');
                ang = ang + ang_blade;
            end

            % Hub(cilinder)
            m=100;
            % create Hub disc
            r_hub = obj.r_bar(1)*obj.R;
            xh=linspace(-r_hub,r_hub,m);
            yh=sqrt((r_hub)^2-xh.^2);
            xh=[xh,flip(xh)]';
            yh=[yh,-yh]';
            Xh=repmat(xh,1,3);
            Yh=repmat(yh,1,3);
            % estrusion of the th Hub Disc
            h_hub  = 0.2*r_hub;
            zh = h_hub*[-1,0, 1 + obj.c(1)*sin(obj.theta(1))/h_hub];
            Zh = repmat(zh,length(xh),1);  
            surf(Xh,Yh,Zh,'FaceColor','k')
            daspect([1 1 1])
            
            xlabel('X [m]')
            ylabel('Y [m]')
%             % Hub(parabolic)
%             m=100;
%             % create Hub disc
%             xh=linspace(-obj.r_bar(1),obj.r_bar(1),m)'*obj.R;
%             yh=xh;
%             [Xh,Yh]=meshgrid(xh,yh);
%             % estrusion of the th Hub Disc
%             k=20;
%             Zh =-k*Xh.^2-k*Yh.^2+0.2;   
%             
%             hold on
%             surf(Xh,Yh,Zh)
%             daspect([1 1 1])
%             camlight HEADLIGHT 
            
        end
        
        
    end
    
end

