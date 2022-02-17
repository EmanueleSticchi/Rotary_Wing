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
        Analisi    
    end
    methods
        %% Auxiliary methods
        % compute ambient conditions
        function obj = ambient(obj)
            [obj.temp, obj.sound_vel, obj.press, obj.rho] = atmosisa(obj.h);
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
        function [P_induced_MR, P_parasite_MR, P_fusolage_MR,...
                P_induced_TR, P_parasite_TR, P_req_MR,...
                P_req_TR, P_req] = Req_power_level_flight(obj)
            %------------------------------------------------------------------
            % This function compute the required power curve in the
            % Power-Velocity plane and stores them inside arrays. The
            % arrays have dimensions 1 x n_vel, the dimension of the
            % velocity vector, that can be changed within the class
            % Input:
            % -
            % Output:
            % - P_induced_MR  : induced power by the main rotor
            % - P_parasite_MR : parasite power of the main rotor
            % - P_fusolage_MR : parasite power of the fusolage and other
            % aerodynamics exposed surfaces
            % - P_induced_TR  : induced power by the tail rotor
            % - P_parasite_TR : parasite power of the tail rotor
            %-------------------------------------------------------------------
            % required thrust
            T_TPP     = obj.W;
            C_T_MR    = T_TPP/( obj.rho*obj.A_D_MR*( obj.omegaR_MR )^2 );
            mu        = obj.V_inf/( obj.omegaR_MR );
            v_i_MR    = C_T_MR./( 2*mu );

            % required power (Main rotor) [W]
            % P_indotta_MR   = W^2./( 2*rho*A_D*V_inf );
            P_induced_MR   = obj.K_MR*T_TPP*v_i_MR/(10^3);
            P_parasite_MR  = obj.Cd_MR*0.125*obj.rho*obj.A_D_MR*obj.V_inf.^3.*( 1 + obj.K1_MR*mu.^2 )/(10^3);
            P_fusolage_MR  = obj.f*0.5*obj.rho*obj.V_inf.^3/(10^3);
            P_req_MR       = P_induced_MR + P_parasite_MR + P_fusolage_MR;

            Q_MR           = P_req_MR*obj.omegaR_TR;
            T_TR           = Q_MR/obj.b;
            CT_TR          = T_TR/( obj.rho*obj.A_D_TR*( obj.omegaR_TR )^2 );
            v_i_TR         = CT_TR./( 2*mu );
            % required power (Tail rotor) [W]
            % P_indotta_TR   = W^2./( 2*rho*A_D*V_inf );
            P_induced_TR   = obj.K_TR*T_TR.*v_i_TR/(10^3);
            P_parasite_TR  = obj.Cd_TR*0.125*obj.rho*obj.A_D_TR*obj.V_inf.^3.*( 1 + obj.K1_TR*mu.^2 )/(10^3);
            P_req_TR       = P_induced_TR + P_parasite_TR;

            % total required power [W]
            P_req = (P_req_MR + P_req_TR + obj.P_req_AUX)*obj.eff_trassm;

        end

    end
end