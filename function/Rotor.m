classdef Rotor
    properties
        % ---------------------------------------------------------------------
        % Geometria
        % ---------------------------------------------------------------------
        N     {mustBeInteger, mustBeFinite}          % Numero di pale
        D     {mustBePositive, mustBeFinite}         % Diametro rotore, [m]
        R     {mustBePositive, mustBeFinite}         % Raggio dell'elica, [m]
        theta (:,1){mustBeReal, mustBeFinite}        % Angoli di calettamento, [rad]
        c     (:,1){mustBeNonnegative, mustBeFinite} % Corda delle sezioni, [m]
        sigma (:,1){mustBeNonnegative, mustBeFinite} % Solidità
        % ---------------------------------------------------------------------
        % Funzionamento
        % ---------------------------------------------------------------------
        RPM   {mustBePositive, mustBeFinite}         % Giri al minuto
        n     {mustBePositive, mustBeFinite}         % Giri al secondo
        omega {mustBePositive, mustBeFinite}         % Velocità di rotazione, [rad/s]
        % ---------------------------------------------------------------------
        % Condizioni ambientali
        % ---------------------------------------------------------------------
        h     {mustBePositive, mustBeFinite} = 0     % Quota di funzionamento
        V_inf (:,1){mustBePositive, mustBeFinite}    % Velocità di funzionamento
        % ---------------------------------------------------------------------
        % Aerodinamica
        % ---------------------------------------------------------------------
        Cl= @(alpha) 2*pi*alpha;
        Cd= @(alpha) 0.01*alpha./alpha;
        % ---------------------------------------------------------------------
        % Analisi
        % ---------------------------------------------------------------------
        Analisi

    end
    properties(SetAccess = private, GetAccess = public)
        % Vettore delle stazioni radiali
        r_bar (:,1) {mustBeInRange(r_bar,0,1,'exclude-lower')} %exclusive
        n_r
        %
        n_analisi = 0
        % condizioni ambientali
        rho   {mustBePositive, mustBeFinite}         % Densità dell'aria
        press {mustBePositive, mustBeFinite}         % pressione dell'aria
        sound_vel {mustBePositive, mustBeFinite}     % velocitò del suono dell'aria
        temp  {mustBePositive, mustBeFinite}         % temperatura dell'aria
    end
    methods
        % imposta il dominio radiale
        function obj = r_(obj,vec_r)
            obj.r_bar=vec_r;
            obj.n_r=length(obj.r_bar);
        end
        % calcolo della solidità
        function obj = sigma_(obj)
            obj.sigma=obj.N/(2*pi)*obj.c.*obj.r_bar.^-1*obj.R;
        end
        % calcolo della condizioni ambientali
        function obj = ambient_(obj,h)
            [obj.temp, obj.sound_vel, obj.press, obj.rho] = atmosisa(h);
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
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Teoria dell'elemento di pala generale per il rotore (SOLVER)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dTc,dQc] = BEMT_rotore(obj,omega,sigma,teta,R,r_bar,V_inf)
            %------------------------------------------------------------------
            % Questa funzione consente di calolcare le prestazioni del rotore
            % una volta fissata la velocità di traslazione, ovvero il parametro
            % mu. La procedura di calcolo assume valida l'ipotei di
            % trascurabilità dell'induzione radiale, che conduce alla
            % definizione di una teoria esplicita per il calcolo delle
            % prestazioni del rotore. La funzione calcola le prestazioni per
            % singolo valore del rapporto di funzionamento, tuttavia può
            % ricevere in input un vettore di rapporti di velocità. In tal caso
            % restituisce una matrice di gradienti di spinta e coppia aventi
            % per dimensioni dimesione_radiale x dimensione vettore velocità
            % Input:
            % - omega: velocità di rotazione del rotore.
            % - sigma: solidità del rotore.
            % - teta:  calettamento del rotore.
            % - R:     raggio del rotore.
            % - r_bar: raggio adimensionalizzato del rotore.
            % - V_inf: velocità di traslazione, necessaria per il calcolo del
            %           rapporto di funzionamento, mu.
            % Output:
            % - dTc/dr: gradiente di spinta per stazione fissata lungo la pala
            % - dQc/dr: gradiente di coppia per stazione fissata lungo la pala
            %-------------------------------------------------------------------
            mu = V_inf./( R*omega );
            B     = mu + ( Cl_alpha.*sigma )/8;
            B2    = B.*B;
            C     = r_bar.*Cl_alpha.*sigma*0.1250.*( teta - (mu./r_bar) );
            lam_i = 0.5*( sqrt( B2 + 4*C ) - B );
            % angolo di inflow
            phi   = ( mu + lam_i )./r_bar;
            % angolo di attacco
            alpha = teta - phi;
            % distribuzione di spinta, coppia
            dTc   = 0.5.*sigma.*obj.Cl(alpha).*(r_bar.^2);
            dQc   = 0.5.*sigma.*(obj.Cl(alpha).*phi + obj.Cd(alpha)).*(r_bar.^3);
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