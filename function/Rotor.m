classdef Rotor
    properties
        % ---------------------------------------------------------------------
        % Geometry
        % ---------------------------------------------------------------------
        N     {mustBeInteger, mustBeFinite}          % Number of blades, [\]
        R     {mustBePositive, mustBeFinite}         % Rotor Radius, [m]
        theta_t (1,1){mustBeReal, mustBeFinite}      % pitch twist angle, [rad]
        c     (:,1){mustBeNonnegative, mustBeFinite} % Rotor chord, [m]
        c_mean{mustBeNonnegative, mustBeFinite}      % Rotor chord, [m]
        sigma {mustBeNonnegative, mustBeFinite}  % Mean solidity, [\]
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
    
    methods(Access=public)
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
            obj = obj.derived_properties();
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
        function obj = BEMT_salita(obj,V_inf,theta0,options)
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
            % - theta0: Comando collettivo,pitch alla radice [rad]
            % - options: Analisys options (see BEMTset_rotor)
            % Output:
            % - dTc/dr: gradiente di spinta per stazione fissata lungo la pala
            % - dQc/dr: gradiente di coppia per stazione fissata lungo la pala
            %-------------------------------------------------------------------
            arguments
                obj
                V_inf      (:,1){mustBeFinite}
                theta0     (1,1){mustBeFinite}
                options = BEMTset_rotor();
            end
            theta = theta0 + obj.theta_t*obj.r_bar;
            for i=1:length(V_inf)
                mu(i)        = V_inf(i)/( obj.R*obj.omega );
                for j=1:obj.n_r   
                    B(i,j)       = mu(i) + ( obj.Cl_alpha.*obj.sigma )/8;
                    B2(i,j)      = B(i,j)*B(i,j);
                    C(i,j)       = obj.r_bar(j)*obj.Cl_alpha*obj.sigma/8*...
                        ( theta(j) - (mu(i)/obj.r_bar(j)) );
                    s.lam_i(i,j) = 0.5*( sqrt( B2(i,j) + 4*C(i,j) ) - B(i,j) );
                    % inflow angle
                    s.phi(i,j)   = ( mu(i) + s.lam_i(i,j) )./obj.r_bar(j);
                    % angle of attack
                    s.alpha(i,j) = theta(j) - s.phi(i,j);
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
            s.mu     = mu;
            s.theta0 = theta0;
            s.theta  = theta;
            obj.n_analisi_salita = obj.n_analisi_salita+1;
            obj.Analisi_salita{obj.n_analisi_salita,1} = s;
        end

        function I = simpsons(obj,f,a,b)
            h=(b-a)/(length(f)-1);
            I= h/3*(f(1)+2*sum(f(3:2:end-2))+4*sum(f(2:2:end-1))+f(end));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% BEMT volo traslato per rotore non rigido (articolato)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = BEMT_articulated(obj,valIN,ToTheta,V_inf_Vec,chi,f,options)
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
            % - valIN           : Flag per scegliere il tipo di risoluzione
            %                     delproblema: "T" -> Spinta fissata;
            %                     "Theta" -> collettivo fissato.
            % - ToTheta[N o rad]: Spinta richiesta o comando collettivo 
            %                     (theta0) a seconda del flag valIN
            % - V_inf[m/s]      : velocità di avanzamento del rotore
            % - Chi  [rad]      : angolo di salita del rotore
            % - f    [visc.area]: prodotto dell'area di riferimento per la
            %                       resistenza dell'intero elicottero.           
            % - options         : Analisys options (see BEMTset_rotor)
            arguments
                obj
                valIN     {mustBeMember(valIN,{'T','Theta',})}
                ToTheta   (1,1) {mustBeFinite}
                V_inf_Vec (:,1){mustBeNonnegative,mustBeFinite}
                chi       {mustBeFinite}
                f         {mustBePositive, mustBeFinite}
                options = BEMTset_rotor();
            end

            switch valIN
                case "T"
                    T  = ToTheta;
                    Tc = T/( obj.rho*(obj.omega*obj.R)^2*obj.A_D );
                case "Theta"
                    theta0 = ToTheta;
                otherwise
                    error('Attezione alla scelta del metodo. Scegliere tra "T" o "Theta"')
            end
            
            
            
            for i = 1:length(V_inf_Vec)
                iter      = 0;        % contatore iterazioni
                iter_cond = 1;        % residuo
                V_inf     = V_inf_Vec(i);
                alpha_TPP = 0;        % valore di primo tentativo
                lam       = 0;        % valore di primo tentativo
                if isequal(valIN,"Theta")
                    mu    = ( V_inf*cos(alpha_TPP) )/( obj.omega*obj.R );
%                     [lam_i,lam,Tc]=initialize_routin(obj,mu,alpha_TPP,theta0,options);
                    Tc = 0.5*(obj.sigma*obj.Cl_alpha)*...
                        (theta0/3*(1+1.5*mu^2) + obj.theta_t/4*(1+mu^2)...
                        - lam/2);
                end
                D_fs      = 0.5*obj.rho*V_inf.^2*f;
                lam_c     = V_inf*sin(chi)/( obj.omega*obj.R );
                
                while iter_cond > options.toll
                    mu    = ( V_inf*cos(alpha_TPP) )/( obj.omega*obj.R );
                                                           
                    if sqrt(mu^2 + lam^2) < 1e-2
                        % nella prima iterazione lam=0 e se mu=0 allora
                        % usiamo la formula dell'induzione in hovering
                        lam_i = sqrt(Tc/2);
                    else
                        lam_i = Tc/( 2*sqrt(mu^2 + lam^2) );
                    end
                    lam      = mu*tan(alpha_TPP) + lam_i;
                    if isequal(valIN,"T")
                        theta0   = ( 3 /( 1 + 1.5*mu.^2 ) )*...
                            ( (2*Tc)/(obj.sigma*obj.Cl_alpha) ...
                            - obj.theta_t/4 - obj.theta_t*( mu.^2 )/4 + 0.5*lam );
                    else
                        Tc = 0.5*(obj.sigma*obj.Cl_alpha)*...
                        (theta0/3*(1+1.5*mu^2) + obj.theta_t/4*(1+mu^2)...
                        - lam/2);
                    end
                    
                    Pc0      = obj.Cd_mean*obj.sigma*( 1 + options.k*mu^2 )/8;
                    Pc       = lam_i*Tc + lam_c*Tc + ...
                        mu*( D_fs/(Tc*obj.rho*obj.omega^2*obj.R^4*pi) )*Tc + Pc0;

                    % flap coeffs.
                    beta0    = obj.gamma*( theta0/8*(1 + mu^2) + ...
                        obj.theta_t/10*(1 + 5*(mu^2)/6) - lam/6 );
                    beta1c   = -2*mu*( (4*theta0/3 + obj.theta_t - lam)...
                        /(1 - 0.5*mu^2) );
                    beta1s   = -4*mu/3*beta0/(1 + 0.5*mu^2);

                    % drag and side force coeffs.
                    % induced drag coeff.
                    Hc_i     = obj.sigma*obj.Cl_alpha*0.5*( theta0*( -beta1c/3 + 0.5*mu*lam ) +...
                        obj.theta_t*( -beta1c/4 + mu*lam/4 ) + 3*lam*beta1c/4 + beta0*beta1s/6 + ...
                        mu*( beta0^2 + beta1c^2 )/4 );
                    % parasite drag coeff.
                    Hc_0     = obj.sigma*obj.Cd_mean*mu/4;
                    % total drag coeff.
                    Hc       = Hc_i + Hc_0;
                    % total lateral force coeff.
                    Yc       = -obj.sigma*obj.Cl_alpha*0.5*...
                        ( theta0*( 3*mu*beta0/4 + beta1s*( 1 + 0.5*3*mu^2 )/3 ) +...
                        obj.theta_t*( 0.5*mu*beta0 + beta1s*( 1 + mu^2 )/4 )...
                        - 3*lam*beta1s/4 + beta0*beta1c*( 1/6 - mu^2 ) - ...
                        0.5*3*mu*lam*beta0 - beta1c*beta1s/4);
                    % modified lam let us compute the variation of lam of the previous
                    % iteration
                    lam_temp = lam;
                    lam      = lam_i + lam_c + mu*( Hc/Tc ) + ...
                        mu*( D_fs/(Tc*obj.rho*obj.omega^2*obj.R^4*pi) );
                    alpha_TPP= atan((lam - Tc/(2*sqrt( mu^2 + lam^2 )))/mu);

                    iter_cond= abs(lam - lam_temp);
                    iter     = iter + 1;
                end
                
                s_art.iter_ART(i)      = iter;
                s_art.iter_cond_ART(i) = iter_cond;
                s_art.lam_Vec(i)       = lam;
                s_art.lam_i(i)         = lam_i;
                s_art.lam_c(i)         = lam_c;
                s_art.mu(i)            = mu;
                s_art.Tc(i)            = Tc;
                s_art.Pc_Vec(i)        = Pc;
                s_art.Hc_i_Vec(i)      = Hc_i;
                s_art.Hc_0_Vec(i)      = Hc_0;
                s_art.Yc_Vec(i)        = Yc;
                s_art.beta0_Vec(i)     = beta0;
                s_art.beta1c_Vec(i)    = beta1c;
                s_art.beta1s_Vec(i)    = beta1s;
                s_art.alpha_TPP_Vec(i) = alpha_TPP;
                s_art.theta0(i)        = theta0;
                s_art.theta(:,i)       = theta0 + obj.theta_t*obj.r_bar;

            end
            s_art.options             = options;
            s_art.alpha_e             = alpha_e(obj,s_art);
            obj.n_analisi_articulated = obj.n_analisi_articulated+1;
            obj.Analisi_articulated{obj.n_analisi_articulated,1} = s_art;

        end
        
        % compute angle of attack for each BE: alpha_e(r_bar,Psi)
        function alpha_e = alpha_e(obj,s)
            Psi=s.options.Psi;
            alpha_e=zeros(obj.n_r,length(Psi),length(s.lam_Vec)); 
            for idxV=1 : length(s.lam_Vec)
                b     =  s.beta0_Vec(idxV) + ...
                         s.beta1c_Vec(idxV)*cos(Psi) +...
                         s.beta1s_Vec(idxV)*sin(Psi);
                b_dot = -s.beta1c_Vec(idxV)*sin(Psi) +...
                         s.beta1s_Vec(idxV)*cos(Psi);
                     
                for i=1:obj.n_r
                    for j=1:length(Psi)
                        alpha_e(i,j,idxV) = s.theta(i,idxV) -...
                                            atan2((s.lam_Vec(idxV)  +...
                                          b_dot(j)*obj.r_bar(i)/obj.omega+...
                                          b(j)*s.mu(idxV)*cos(Psi(j))),(...
                                          obj.r_bar(i) + s.mu(idxV)*sin(Psi(j))));
%                         if abs(alpha_e(i,j,idxV)) > 1e2*pi/180
%                             alpha_e(i,j,idxV) = s.theta(i,idxV) - pi/2;
%                         end
                    end
                end
            end
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
        function alphamap(obj,valIN,val)
            % Plot alpha_e contour
            % INPUT:
            % - valIN:    flag per l'input val:
            %                   - 'Plot'  -> in tal caso val dovrà essere
            %                   un cell array 2x1 in cui vi è 
            %                   una struct (output di BEMT_articulated) ed
            %                   un vettore di valori di mu per i quali si
            %                   desidera effettuare i plot
            %                    
            %                   - 'Solve' -> in tal caso val dovrà essere
            %                   una cell array  6x1 con gli input da dare
            %                   alla funzione BEMT_articulated 
            %                 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            switch valIN
                case 'Plot'
                    s=val{1,1};
                    for i =1:length(val{2,1})
                        [~,idxMu(i)]=min(abs(s.mu - val{2,1}(i)));
                    end
                case 'Solve'
                    valIN     = val{1,1};
                    ToTheta   = val{2};
                    V_inf_Vec = val{3};
                    chi       = val{4};
                    f         = val{5};
                    options   = val{6};
                    obj2=obj.BEMT_articulated(valIN,ToTheta,V_inf_Vec,...
                        chi,f,options);
                    s=obj2.Analisi_articulated{obj2.n_analisi_articulated,1};
                    idxMu=1:length(s.mu);
                otherwise
                    error('Attenzione valIN può essere: "Plot" o "Solve"')
            end
            

            for i = 1:length(idxMu)
                figure
                idxV=idxMu(i);
                alpha_e=s.alpha_e(:,:,idxV)*180/pi;
                % Create polar data
                [r,psi] = meshgrid(obj.r_bar,s.options.Psi);
                % Convert to Cartesian
                x = r.*cos(psi);
                y = r.*sin(psi);
                % define polar axes
                h = polar(x,y);
                hold on;
                polar(s.options.Psi,obj.r_bar(1)*ones(length(s.options.Psi),1),'k')
                polar(s.options.Psi,obj.r_bar(end)*ones(length(s.options.Psi),1),'k')
                % contourf(x,y,alpha_e');
                pc= pcolor(x,y,alpha_e');
                contour(x,y,alpha_e','k','ShowText','on');
                shading interp
                % colormap 'hsv'
                cbar=colorbar(gca);
                cbar.Label.String = '\alpha_e';
                cbar.Label.FontSize= 16;
                % cbar.Limits = [-10 10];
  
                % Hide the POLAR function data and leave annotations
                set(h,'Visible','off')
                % Turn off axes and set square aspect ratio
                axis off
                axis image
                view([90 90])
                title(['\mu = ',num2str(s.mu(idxV)),'   \alpha_{e_{max}} = ',...
                    num2str(max(alpha_e,[],'all')),' deg'])
            end
            




        end
        

        
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
    methods(Access=private)
        function [lam_i,lam,Tc]=initialize_routin(obj,mu,alpha_TPP,theta0,options)
            lam=0;
            res=1;
            while res > options.toll
                lam_old=lam;
                Tc = max([0.5*(obj.sigma*obj.Cl_alpha)*...
                    (theta0/3*(1+1.5*mu^2) + obj.theta_t/4*(1+mu^2)...
                    - lam/2),0]);
                lam_i = Tc/( 2*sqrt(mu^2 + lam^2) );
                if isinf(lam_i)
                    % nella prima iterazione lam=0 e se mu=0 allora
                    % usiamo la formula dell'induzione in hovering
                    lam_i = sqrt(Tc/2);
                end
                lam = lam_i +mu*tan(alpha_TPP);
                res = abs(lam -lam_old);
            end
        end
        
    end
end