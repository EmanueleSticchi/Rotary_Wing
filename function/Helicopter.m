classdef Helicopter<Rotor
    % TODO:
    % -> INSERT LIMITS FOR h, fuel_load (0,1) and others    when possible
    properties
        % ---------------------------------------------------------------------
        % Flight Conditions
        % ---------------------------------------------------------------------
        h        {mustBeFinite}                     % [m]
        % cancellare se non la usiamo
        V_max    {mustBePositive, mustBeFinite}     % [m/s]
        % ---------------------------------------------------------------------
        % Mass
        % ---------------------------------------------------------------------
        % poi vediamo
        fuel_load {mustBePositive, mustBeFinite, mustBeInRange(fuel_load,0,1,'exclude-lower')}
        W_fuel    {mustBePositive, mustBeFinite}    % [N]
        W_mtow    {mustBePositive, mustBeFinite}    % [N]
        % ---------------------------------------------------------------------
        % Geometry
        % ---------------------------------------------------------------------
        % Main Rotor
        MR = Rotor();
        % Tail rotor
        TR = Rotor();
        % distance between rotor axes
        lr         {mustBePositive, mustBeFinite}       % [m]
        % ---------------------------------------------------------------------
        % Propulsion
        % ---------------------------------------------------------------------
        engine_number {mustBePositive, mustBeFinite}  % [\]
        SFC           {mustBePositive, mustBeFinite}  % [Kg/kW*h]?
        engine_power  {mustBePositive, mustBeFinite}  % [W]
        % ---------------------------------------------------------------------
        % Power
        % ---------------------------------------------------------------------
        P_req_AUX    {mustBePositive, mustBeFinite}   % [W]
        % loss trasmission coefficient
        eta_t        {mustBePositive, mustBeFinite,...
            mustBeGreaterThan(eta_t,1)}               % [\]
        % ---------------------------------------------------------------------
        % Correction coeffs
        % ---------------------------------------------------------------------
        % fattore di correzione dovuto alla non-uniformita' dell'induzione
        % sul rotore reale: Pc_i = k*lam_i*T_c
        k_i_MR  = 1.2
        k_i_TR  = 1.4
        % Fattore di scorrimento per il calcolo della potenza parassita
        % Pc0 = sigma*Cd_mean/8*(1+k*mu^2)
        k_mu_MR = 4.7
        k_mu_TR = 4.7
    end
    properties(SetAccess = private, GetAccess = public)
        % Access denied: there's no way fo the user to change the value of
        % these variables.

        % Ambient conditions
        rho       {mustBePositive, mustBeFinite}% Air density
        press     {mustBePositive, mustBeFinite}% pressione dell'aria
        sound_vel {mustBePositive, mustBeFinite}% velocità del suono dell'aria
        temp      {mustBePositive, mustBeFinite}% temperatura dell'aria
        n_analisi = 0;
        % Mass
        M         {mustBePositive, mustBeFinite}
        W         {mustBePositive, mustBeFinite}       % [N]
        % ---------------------------------------------------------------------
        % Analisi
        % ---------------------------------------------------------------------
        % Power Analysis
        PA  
        % Number of power analisys
        n_PA
    end
    methods
        %% Auxiliary methods
        % compute ambient conditions
        function obj = ambient(obj,h)
            obj.h = h;
            [obj.temp, obj.sound_vel, obj.press, obj.rho] = atmosisa(obj.h);
%             T0  = 288.15;
%             mu0 = 1.79e-5;
%             obj.mu_visc=mu0*(obj.temp/T0)^1.5*((T0+110)/(obj.temp+110));
        end
        %% Solver methods
        function obj = Required_Power(obj,h,V_inf_Vec,Chi,T,f)
            arguments
                obj
                h         {mustBeNonnegative,mustBeFinite}
                V_inf_Vec (:,1){mustBeNonnegative,mustBeFinite}
                Chi       {mustBeFinite}
                T         {mustBePositive,mustBeFinite}
                f         {mustBePositive, mustBeFinite}
            end
            V_inf_Vec = sort(V_inf_Vec); flag = sum(V_inf_Vec ~= 0);
            V_inf_Vec = V_inf_Vec(V_inf_Vec ~= 0);
            %--------------------------------------------------------------
            % Main Rotor
            obj.MR.h     = h;                   % set altitude
            obj.MR       = obj.MR.ambient();    % compute ambient properties
            options      = BEMTset_rotor();
            options.k_i  = obj.k_i_MR;
            options.k_mu = obj.k_mu_MR;
            obj.MR       = obj.MR.BEMT_articulated('T',T,V_inf_Vec,Chi,f,options);
            s            = obj.MR.Analisi_articulated{obj.MR.n_analisi_articulated,1};
            P_MR         = s.Pc_vec * obj.rho*pi*obj.MR.R^5*obj.MR.omega^3;
            Pi_MR        = s.Pci_vec * obj.rho*pi*obj.MR.R^5*obj.MR.omega^3;
            P0_MR        = s.Pc0_vec * obj.rho*pi*obj.MR.R^5*obj.MR.omega^3;
            P_fus        = s.Pc_fus_Vec * obj.rho*pi*obj.MR.R^5*obj.MR.omega^3;
            Q_MR         = P_MR/obj.MR.omega;
            % Tail Rotor
            T_TR         = Q_MR/obj.lr;
            Tc_TR        = T_TR/( obj.rho*pi*obj.TR.R^4*obj.MR.omega^2 );
            v_i_TR       = sqrt( -0.5*(V_inf_Vec.^2 + sqrt(V_inf_Vec.^4 + ...
                            4*(T_TR.*0.5/obj.rho/(pi*obj.TR.R^2))^4)));
            Pi_TR        = obj.ki_TR*T_TR*v_i_TR;
            mu_TR        = V_inf_vec/obj.TR.omega/obj.TR.R;
            P0_TR        = (obj.TR.Cd_mean*obj.TR.sigma*(1 + obj.k_mu_TR*mu_TR.^2)/8)...
                            *(obj.rho*pi*obj.TR.R^5*obj.TR.omega^3);
            P_TR         = Pi_TR + P0_TR;
            if flag ~= 0
                % nel caso in cui venga richiesta una condizione di hover
                % To do : MR in hover e controllare BEMT_Artic.. in FF
                % livellato
            end
        end
        
        
        
        
        
        
        % required power for level flight
        function obj = Req_power_level_flight(obj,h,V_inf_Vec,Chi,T,f,valIN,PoVc)
            %------------------------------------------------------------------
            % This function compute the required power curve in the
            % Power-Velocity plane and stores them inside arrays. The
            % arrays have dimensions 1 x n_vel, the dimension of the
            % velocity vector, that can be changed within the class
            % Input:
            % - h
            % - V_inf_Vec
            % - Chi
            % Output:
            % - P_induced_MR  : induced power by the main rotor
            % - P_parasite_MR : parasite power of the main rotor
            % - P_fusolage_MR : parasite power of the fusolage and other
            % aerodynamics exposed surfaces
            % - P_induced_TR  : induced power by the tail rotor
            % - P_parasite_TR : parasite power of the tail rotor
            %-------------------------------------------------------------------
            arguments
                obj
                h         {mustBeNonnegative,mustBeFinite}
                V_inf_Vec (:,1){mustBeNonnegative,mustBeFinite}
                Chi       {mustBeFinite}
                T         {mustBePositive,mustBeFinite}
                f         {mustBePositive, mustBeFinite}
                valIN     {mustBeMember(valIN,{'P','Chi'})}
                PoVc     (1,1) {mustBeFinite,mustBeMember}
            end
            
            % required thrust
            T_TPP     = T;
            Tc_MR    = T_TPP/( obj.rho*pi*obj.MR.R^4*obj.MR.omega^2 );
            mu       = V_inf_Vec/( obj.MR.omega*obj.MR.R );
            lam_i_MR = sqrt( -0.5*(V_inf_Vec.^2 + sqrt(V_inf_Vec.^4 + ...
                            4*(T_TPP.*0.5/obj.rho/(pi*obj.MR.R^2))^4)));
            D_fs     = 0.5*obj.rho*V_inf_Vec.^2*f;

            
            % required power (Main rotor) [W]
            s.Pc_induced_MR   = obj.k_i_MR*T_TPP*lam_i_MR;
            s.Pc_parasite_MR  = obj.MR.sigma*obj.MR.Cd_mean*( 1 + obj.k_mu_MR*mu.^2 )/8;
            s.Pc_fusolage_MR  = mu*( D_fs/T_TPP )*Tc_MR;
            s.Pc_req_MR       = s.Pc_induced_MR + s.Pc_parasite_MR + s.Pc_fusolage_MR;


            Q_MR           = s.P_req_MR/obj.MR.omega;
            T_TR           = Q_MR/obj.b;
            Tc_TR          = T_TR/( obj.rho*pi*obj.MR.R^4*obj.MR.omega^2 );
            lam_i_TR       = sqrt( -0.5*(V_inf_Vec.^2 + sqrt(V_inf_Vec.^4 + ...
                            4*(T_TR.*0.5/obj.rho/(pi*obj.TR.R^2))^4)));
            % required power (Tail rotor) [W]
            s.Pc_induced_TR   = obj.k_i_TR*T_TR*lam_i_TR;
            s.Pc_parasite_TR  = obj.TR.sigma*obj.TR.Cd_mean*( 1 + obj.k_mu_TR*mu.^2 )/8;
            s.Pc_req_TR       = s.P_induced_TR + s.P_parasite_TR;

            % total required power [W]       
            s.Pi_MR        = s.Pc_induced_MR * obj.rho*pi*obj.MR.R^5*obj.MR.omega^3;
            s.P0_MR        = s.Pc_parasite_MR * obj.rho*pi*obj.MR.R^5*obj.MR.omega^3;
            s.P_fus_MR     = s.Pc_fusolage_MR * obj.rho*pi*obj.MR.R^5*obj.MR.omega^3;
            s.P_req_MR     = s.Pc_req_MR * obj.rho*pi*obj.MR.R^5*obj.MR.omega^3;

            s.Pi_TR        = s.Pc_induced_TR * obj.rho*pi*obj.TR.R^5*obj.TR.omega^3;
            s.P0_TR        = s.Pc_parasite_TR * obj.rho*pi*obj.TR.R^5*obj.TR.omega^3;
            s.P_req_TR     = s.Pc_req_TR * obj.rho*pi*obj.TR.R^5*obj.TR.omega^3;

            s.P_req_hori   = (s.P_req_MR + s.P_req_TR + obj.P_req_AUX)*obj.eta_t;

            switch valIN
                case 'P' % available power
                    Vc         = (PoVc - s.P_req_hori)/TPP;
                    s.Pc_climb = Vc*TPP;
                case 'Vc'
                    Vc         = PoVc;
                    s.Pc_climb = Vc*TPP;
            end
            
            obj.n_PA = obj.n_PA+1;
            obj.PA{obj.n_PA,1} = s;

        end

    end
end