classdef Rotor
    properties
        % ---------------------------------------------------------------------
        % Geometry
        % ---------------------------------------------------------------------
        N     {mustBeInteger, mustBeFinite}          % Number of blades, [\]
        R     {mustBePositive, mustBeFinite}         % Rotor Radius, [m]
        theta (:,1){mustBeReal, mustBeFinite}        % pitch setting, [rad]
        c     (:,1){mustBeNonnegative, mustBeFinite} % Rotor chord, [m]
        c_mean{mustBeNonnegative, mustBeFinite}      % Rotor chord, [m]
        sigma{mustBeNonnegative, mustBeFinite}  % Mean solidity, [\]
        I     {mustBePositive, mustBeFinite}         % Moment of inertia, [kgm^2]
        gamma {mustBePositive, mustBeFinite}         % Lock number
        % ---------------------------------------------------------------------
        % Working conditions
        % ---------------------------------------------------------------------
        RPM   {mustBePositive, mustBeFinite}         % Giri al minuto
        n     {mustBePositive, mustBeFinite}         % Giri al secondo
        omega {mustBePositive, mustBeFinite}         % Velocità di rotazione, [rad/s]
        Mach_limit = 0.7                             % Maximum Mach number at the tip of the blade
        % ---------------------------------------------------------------------
        % Ambient conditions
        % ---------------------------------------------------------------------
        h     {mustBeNonnegative, mustBeFinite}       % Quota di funzionamento
        % ---------------------------------------------------------------------
        % Aerodynamics
        % ---------------------------------------------------------------------
        Cl = @(alpha) 2*pi*alpha;
        Cd = @(alpha) 0.01*alpha./alpha;
        Cl_alpha = 2*pi;
        Cd_mean  = 0.01;
        % ---------------------------------------------------------------------
        % Analisys & Design
        % ---------------------------------------------------------------------


    end
    properties(SetAccess = private, GetAccess = public)
        % radial vector
        r_bar (:,1) {mustBeInRange(r_bar,0,1,'exclude-lower')} %exclusive
        n_r
        %
        n_analisi_salita = 0
        n_analisi_articulated = 0
        n_analisi_autorot = 0
        % ambient conditons
        rho       {mustBePositive, mustBeFinite}
        press     {mustBePositive, mustBeFinite}
        sound_vel {mustBePositive, mustBeFinite}
        temp      {mustBePositive, mustBeFinite}
        % geometry
        D     {mustBePositive, mustBeFinite}         % Rotor diameter, [m]
        A_D   {mustBePositive, mustBeFinite}         % Rotor area, [m^2]

        % storage variables for analysis and design
        Analisi_salita
        Analisi_articulated;
        Analisi_autorot;
        Design
    end
    methods
        % Compute ambient conditions. This function needs to be called
        % first in order to set atmospheric conditions properly.
        function obj = ambient(obj)
            [obj.temp, obj.sound_vel, obj.press, obj.rho] = atmosisa(obj.h);
        end

        % set radial domain
        function obj = r(obj,vec_r)
            obj.r_bar = vec_r;
            obj.n_r   = length(obj.r_bar);
        end
        
        % mass properties
        function obj = mass_prop(obj,valIN,val)
            obj = obj.derived_properties;
            switch valIN
                case 'G'
                    obj.gamma = val;
                    obj.I     = obj.Cl_alpha*obj.rho*obj.R^4*obj.c_mean/val;
                case 'I'
                    obj.I= val;
                    obj.gamma = obj.Cl_alpha*obj.rho*obj.R^4*obj.c_mean./obj.I;
                otherwise
                    mustBeMember(valIN,{'G','I'})
            end
        end
        
        % Compute some mass and geometric property. This function needs to
        % be called before doing any other calculations but still after the
        % definition of the main properties. The function needs to be
        % called only once. The unction computes any derived property.
        function obj = derived_properties(obj)
            obj.D     = obj.R*2;
            obj.A_D   = pi*obj.R*obj.R;
            obj.c_mean= mean(obj.c);
            obj.sigma = ( obj.c_mean*obj.N )/( pi*obj.R );
        end

        % Compute rotational velocity. This function let the user decide
        % the tip speed depending on the input parameters: sometime the
        % rotational velocity of the rotor can be known, otherwise the
        % upper limit for the tip speed is set by compressibility
        % requirements: in this case the user should select the value 1 for
        % the flag. While the flag is equal to 1 there's no need to give
        % the function other inputs.
        function obj = rot_vel(obj,valIN,val,flag_external_tip_speed)
            arguments
                obj
                valIN
                val
                flag_external_tip_speed = 0
            end
            if flag_external_tip_speed == 1
                obj.omega = obj.sound_vel*obj.Mach_limit/obj.R;
                obj.n     = obj.omega/(2*pi);
                obj.RPM   = obj.n*60;
            else
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
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% BEMT salita assiale
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = BEMT_salita(obj,V_inf,options)
            %------------------------------------------------------------------
            % Questa funzione consente di calcolare le prestazioni del rotore
            % una volta fissata la velocità di salita, ovvero il parametro
            % mu. La procedura di calcolo assume valida l'ipotesi di
            % trascurabilità dell'induzione radiale, che conduce alla
            % definizione di una teoria esplicita per il calcolo delle
            % prestazioni del rotore. La funzione calcola le prestazioni per
            % singolo valore del rapporto di funzionamento, tuttavia può
            % ricevere in input un vettore di rapporti di velocità. In tal caso
            % restituisce una matrice di gradienti di spinta e coppia aventi
            % per dimensioni dimesione_radiale x dimensione vettore velocità
            % Input:
            % - obj: 
            %       - omega: velocità di rotazione del rotore.
            %       - sigma: solidità del rotore.
            %       - teta:  calettamento del rotore.
            %       - R:     raggio del rotore.
            %       - r_bar: raggio adimensionalizzato del rotore.
            % - V_inf: velocità di traslazione, necessaria per il calcolo del
            %           rapporto di funzionamento, mu.
            % - options: Analisys options (see BEMTset_rotor)
            % Output:
            % - dTc/dr: gradiente di spinta per stazione fissata lungo la pala
            % - dQc/dr: gradiente di coppia per stazione fissata lungo la pala
            %-------------------------------------------------------------------
            arguments
                obj
                V_inf      (:,1){mustBeFinite}
                options = BEMTset_rotor();
            end
            for i=1:length(V_inf)
                mu(i)        = V_inf(i)/( obj.R*obj.omega );
                for j=1:obj.n_r   
                    B(i,j)       = mu(i) + ( obj.Cl_alpha.*obj.sigma )/8;
                    B2(i,j)      = B(i,j)*B(i,j);
                    C(i,j)       = obj.r_bar(j)*obj.Cl_alpha*obj.sigma/8*( obj.theta(j) - (mu(i)/obj.r_bar(j)) );
                    s.lam_i(i,j) = 0.5*( sqrt( B2(i,j) + 4*C(i,j) ) - B(i,j) );
                    % inflow angle
                    s.phi(i,j)   = ( mu(i) + s.lam_i(i,j) )./obj.r_bar(j);
                    % angle of attack
                    s.alpha(i,j) = obj.theta(j) - s.phi(i,j);
                    s.Cl(i,j)    = obj.Cl(s.alpha(i,j));
                    s.Cd(i,j)    = obj.Cd(s.alpha(i,j));
                    % thrust and torque distributions
                    s.dTc(i,j)   = 0.5*obj.sigma*s.Cl(i,j)*(obj.r_bar(j)^2);
                    s.dQc(i,j)   = 0.5*obj.sigma*(s.Cl(i,j)*s.phi(i,j) + s.Cd(i,j))*(obj.r_bar(j)^3);
                end

                % thrust and torque
                s.Tc(i) = obj.simpsons(s.dTc(i,:),obj.r_bar(1),options.B*obj.r_bar(end));
                s.Qc(i) = obj.simpsons(s.dQc(i,:),obj.r_bar(1),obj.r_bar(end));

            end
            s.mu = mu;
            obj.n_analisi_salita = obj.n_analisi_salita+1;
            obj.Analisi_salita{obj.n_analisi_salita,1} = s;
        end

        function I = simpsons(obj,f,a,b)
            h=(b-a)/(length(f)-1);
            I= h/3*(f(1)+2*sum(f(3:2:end-2))+4*sum(f(2:2:end-1))+f(end));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Analisi del rotore articolato (coefficienti di flappeggio)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = articulated_rotor(obj,V_inf_Vec,chi,f,W,theta_t)
            %------------------------------------------------------------------
            % Questa funzione consente di calcolare i coefficienti adimensionali di spinta,
            % resistenza, forza laterale, coppia e potenza, gli angoli di
            % flappeggio, l'angolo di attacco del rotore e l'angolo di
            % inflow. Tutte le operazioni si basa sulle ipotesi di piccoli
            % angoli, sono richiesti: piccoli angoli di attacco del rotore,
            % piccoli angoli di flappeggio, piccoli angoli di salita.
            % Vengono poi supposte costanti le velocità di traslazione dell'elicottero e
            % angolare del rotore. Si assume induzione costante sul rotore
            % e si considera come unica cerniera presente quella di
            % flappeggio (cfr. ipotesi per il rotore articolato). Questa
            % funzione è stata progettata anche per lavorare in coppia con
            % la classe Helicopter.m, tuttavia qualora fornite le variabili
            % di input è possibile utilizzare questo metodo come standalone
            % in un codice "ad hoc", in cui si utiizza la sola classe
            % Rotor.m.
            % Input:
            % - V_inf[m/s]      : velocità di avanzamento del rotore
            % - Chi  [rad]      : angolo di salita del rotore
            % - f    [visc.area]: prodotto dell'area di riferimento per la
            % resistenza dell'intero elicottero.
            % - W    [N]        : peso dell'elicottero
            % - theta_w [rad]   : coefficiente angolare della retta che
            % definisce lo svergolamento
            arguments
                obj
                V_inf_Vec (:,1){mustBeNonnegative,mustBeFinite}
                chi       {mustBeFinite}
                f         {mustBePositive, mustBeFinite}
                W         {mustBePositive, mustBeFinite}
                theta_t   {mustBeFinite}
            end
            if chi > 0.1745
                warning('The input value for the climb angle may not be small enough!')
            end

            alpha_TPP     = 0;
            toll          = 1e-5;
            iter          = 0;
            s_art.Tc      = W/( obj.rho*(obj.omega*obj.R)^2*obj.A_D );
            for i = 1:length(V_inf_Vec)
                V_inf     = V_inf_Vec(i);
                D_fs      = 0.5*V_inf.^2*f;
                lam_c     = V_inf*sin(chi)/( obj.omega*obj.R );
                iter_cond = 1;
                while iter_cond > toll

                    mu = ( V_inf*cos(alpha_TPP) )/( obj.omega*obj.R );
                    if mu <= 0.1    % mu < 0.1 || mu = 0.1
                        if iter == 0% during the first iteration lam it's not defined,
                            % we use an approximated value for the induction. In this case
                            % the induction is computed with the hovering formula since mu
                            % < 0.1. It is important to notice that this value it's used
                            % only for the first iteration.
                            %               lam_i = sqrt(W/2*rho*A_D)/OmegaR;
                            lam_i = sqrt(-V_inf^2/2 + sqrt( V_inf^4/4 + (W/(2*obj.rho*obj.A_D))^2 ))/(obj.omega*obj.R);
                        else
                            lam_i = s_art.Tc/( 2*sqrt(mu.^2 + lam.^2) );
                        end
                    end

                    if mu > 0.1     % mu > 0.1
                        if iter == 0% during the first iteration lam it's not defined,
                            % we use an approximated value for the induction.
                            lam_i = s_art.Tc/( 2*mu );
                        else
                            lam_i = s_art.Tc/( 2*sqrt(mu.^2 + lam.^2) );
                        end
                    end

                    lam      = mu*tan(alpha_TPP) + lam_i;
                    theta0   = ( 3 /( 1 + 1.5*mu.^2 ) )*( (2*s_art.Tc)/(obj.sigma*obj.Cl_alpha) ...
                        - theta_t/4 - theta_t*( mu.^2 )/4 + 0.5*lam );

                    Pc0      = obj.Cd_mean*obj.sigma*( 1 + 3*mu.^2 )/8;
                    Pc       = lam_i*s_art.Tc + lam_c*s_art.Tc + mu*( D_fs/W )*s_art.Tc + Pc0;

                    % flap coeffs.
                    beta0    = obj.gamma*( theta0/8*(1 + mu.^2) + ...
                        theta_t/10*(1 + 5*(mu.^2)/6) - lam/6 );
                    beta1c   = -2*mu*( (4*theta0/3 + theta_t - lam)...
                        /(1 - 0.5*mu.^2) );
                    beta1s   = -4*mu/3.*beta0/(1 + 0.5*mu.^2);

                    % drag and side force coeffs.
                    % induced drag coeff.
                    Hc_i     = obj.sigma*obj.Cl_alpha*0.5*( theta0*( -beta1c/3 + 0.5*mu*lam ) +...
                        theta_t*( -beta1c/4 + mu*lam/4 ) + 3*lam*beta1c/4 + beta0*beta1s/6 + ...
                        mu*( beta0^2 + beta1c^2 )/4 );
                    % parasite drag coeff.
                    Hc_0     = obj.sigma*obj.Cd_mean*mu/4;
                    % total drag coeff.
                    Hc       = Hc_i + Hc_0;
                    % total lateral force coeff.
                    Yc       = -obj.sigma*obj.Cl_alpha*0.5*( theta0*( 3*mu*beta0/4 + beta1s*( 1 + 0.5*3*mu^2 )/3 ) +...
                        theta_t*( 0.5*mu*beta0 + beta1s*( 1 + mu^2 )/4 ) - 3*lam*beta1s/4 + beta0*beta1c*( 1/6 - mu^2 ) - ...
                        0.5*3*mu*lam*beta0 - beta1c*beta1s/4);
                    % modified lam let us compute the variation of lam of the previous
                    % iteration
                    lam_temp = lam;
                    lam      = lam_i + lam_c + mu*( Hc/s_art.Tc ) + mu*( D_fs/W );
                    alpha_TPP= atan((lam - s_art.Tc/(2*sqrt( mu^2 + lam^2 )))/mu);

                    iter_cond= abs(lam - lam_temp);
                    iter     = iter + 1;
                end
                s_art.iter_ART(i)      = iter;
                s_art.iter_cond_ART(i) = iter_cond;
                s_art.lam_Vec(i)       = lam;
                s_art.Pc_Vec(i)        = Pc;
                s_art.Hc_i_Vec(i)      = Hc_i;
                s_art.Hc_0_Vec(i)      = Hc_0;
                s_art.Yc_Vec(i)        = Yc;
                s_art.beta0_Vec(i)     = beta0;
                s_art.beta1c_Vec(i)    = beta1c;
                s_art.beta1s_Vec(i)    = beta1s;
                s_art.alpha_TPP_Vec(i) = alpha_TPP;
                if alpha_TPP > 0.1745 % this warning is inserted just to highlight the fact that
                    % the angle of attack of the rotor it's greater than 10°.
                    warning('The rotor angle of attack may not be small enough!')
                end

            end
            obj.n_analisi_articulated = obj.n_analisi_articulated+1;
            obj.Analisi_articulated{obj.n_analisi_articulated,1} = s_art;

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Analisi dell'autorotazione 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = autorot_rotor(obj,V_inf_Vec,chi,f,W,theta_t)
            %------------------------------------------------------------------
            % Questa funzione consente di calcolare le caratteristiche
            % della manovra di autorotazione per il rotore in analisi. Come
            % descritto a pagina 110 degli appunti, è necessario imporre i
            % parametri lambda_climb e mu. Come prescritto dalla manovra di
            % autorotazione si impone che P = 0 e si ricava una relazione
            % per theta_0. Noto lambda_climb si possono calcolare lambda e
            % lambda indotto. Noti quesri coefficienti si possono calcolare
            % i coefficienti di forza. Noto Tc e il peso W si può ricavare
            % Omega, la velocità angolare del rotore. Infine si calcola
            % alfa e V_infinito.

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PLOTTING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Model3D(obj,x,z)
            % Plot 3D model of the propeller
            % INPUT:
            % - x e z               Airfoil coordinate (in percent of chord)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                x_rot=data_rot(:,1)*obj.c(i);
                z_rot=data_rot(:,2)*obj.c(i);
                X(:,i)=x_rot-mean(x_rot);
                Z(:,i)=z_rot;
            end

            Y=obj.r_bar'.*obj.R;
            Y=repmat(Y,length(x),1);
            light('Style','local','Position',[1 -1 0]);

            s=surf(X,Y,Z,'FaceColor',[0.65 0.65 0.65],'FaceLighting',...
                'gouraud','EdgeColor','none');
            daspect([1 1 1])
            camlight right
            % material dull

            % Hub(cilinder)
            m=100;
            % create Hub disc
            xh=linspace(-obj.r_bar(1),obj.r_bar(1),m)*obj.R;
            yh=sqrt((obj.r_bar(1)*obj.R)^2-xh.^2);
            xh=[xh,flip(xh)]';
            yh=[yh,-yh]';
            Xh=repmat(xh,1,3);
            Yh=repmat(yh,1,3);
            % estrusion of the th Hub Disc
            Zh=repmat(-0.1:0.1:0.1,length(xh),1);
            hold on
            surf(Xh,Yh,Zh)
            daspect([1 1 1])

            % Hub(parabolic)
            m=100;
            % create Hub disc
            xh=linspace(-obj.r_bar(1),obj.r_bar(1),m)'*obj.R;
            yh=xh;
            [Xh,Yh]=meshgrid(xh,yh);
            % estrusion of the th Hub Disc
            k=20;
            Zh=-k*Xh.^2-k*Yh.^2+0.2;
            hold on
            surf(Xh,Yh,Zh)
            daspect([1 1 1])
            camlight HEADLIGHT
            % other blade
            ang_blade=2*pi/obj.N;
            ang=0;

            for i =2:obj.N
                ang=ang+ang_blade;
                xyz=[X(:),Y(:),Z(:)];
                xyz_rot=xyz*R3(ang);
                X_rot=reshape(xyz_rot(:,1),[length(x),obj.n_r]);
                Y_rot=reshape(xyz_rot(:,2),[length(x),obj.n_r]);
                Z_rot=reshape(xyz_rot(:,3),[length(x),obj.n_r]);
                s=surf(X_rot,Y_rot,Z_rot,'FaceColor',[0.65 0.65 0.65],...
                    'FaceLighting','gouraud','EdgeColor','none');

            end
        end
    end
end