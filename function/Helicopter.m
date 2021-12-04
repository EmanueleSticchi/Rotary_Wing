classdef Helicopter<Rotor
    % TODO:
    % -> INSERT LIMITS FOR h, fuel_load (0,1) and otherswhen possible
    properties
        % ---------------------------------------------------------------------
        % Flight Conditions
        % ---------------------------------------------------------------------
        h        {mustBeFinite}   % [m]
        V_max    {mustBePositive, mustBeFinite}   % [m/s]
        % ---------------------------------------------------------------------
        % Mass
        % ---------------------------------------------------------------------
        fuel_load {mustBePositive, mustBeFinite}
        M_fuel    {mustBePositive, mustBeFinite}    % [kg]
        M_nofuel  {mustBePositive, mustBeFinite}    % [kg]
        % ---------------------------------------------------------------------
        % Geometry
        % ---------------------------------------------------------------------
        % Main Rotor
        %     Main_rotor = Rotor();
        R_MR      {mustBePositive, mustBeFinite}  % [m]
        c_MR      {mustBePositive, mustBeFinite}  % [m]
        N_MR      {mustBePositive, mustBeFinite}  % [\]
        sigma_MR  {mustBePositive, mustBeFinite}  % [\]
        % Tail rotor
        %     Tail_rotor = Rotor();
        R_TR      {mustBePositive, mustBeFinite}  % [m]
        c_TR      {mustBePositive, mustBeFinite}  % [m]
        N_TR      {mustBePositive, mustBeFinite}  % [\]
        sigma_TR  {mustBePositive, mustBeFinite}  % [\]
        % ---------------------------------------------------------------------
        % Aerodynamics
        % ---------------------------------------------------------------------
        % Main Rotor
        RPM_MR   {mustBePositive, mustBeFinite}   % [rpm]
        n_MR     {mustBePositive, mustBeFinite}   % [1/s]
        omega_MR {mustBePositive, mustBeFinite}   % [rad/s]
        % Tail rotor
        RPM_TR   {mustBePositive, mustBeFinite}   % [rpm]
        n_TR     {mustBePositive, mustBeFinite}   % [1/s]
        omega_TR {mustBePositive, mustBeFinite}   % [rad/s]
        Cl = @(alpha) 2*pi*alpha;
        Cd = @(alpha) 0.01*alpha./alpha;
        K1_MR    {mustBePositive, mustBeFinite}
        K_MR     {mustBePositive, mustBeFinite}
        f        {mustBePositive, mustBeFinite}
        K1_TR    {mustBePositive, mustBeFinite}
        K_TR     {mustBePositive, mustBeFinite}
        Cd_MR    {mustBePositive, mustBeFinite}
        Cd_TR    {mustBePositive, mustBeFinite}
        % ---------------------------------------------------------------------
        % Propulsion
        % ---------------------------------------------------------------------
        engine_number {mustBePositive, mustBeFinite}  % [\]
        engine_power  {mustBePositive, mustBeFinite}  % [W]
        % ---------------------------------------------------------------------
        % Power
        % ---------------------------------------------------------------------
        P_req_AUX    {mustBePositive, mustBeFinite}
        eff_trassm   {mustBePositive, mustBeFinite}
        % ---------------------------------------------------------------------
        % Analisi
        % ---------------------------------------------------------------------
        Analisi

    end
    properties(SetAccess = private, GetAccess = public)
        % Access denied: there's no way fo the user to change the value of
        % these variables.

        % Ambient conditions
        rho       {mustBePositive, mustBeFinite}% Air density
        press     {mustBePositive, mustBeFinite}% pressione dell'aria
        sound_vel {mustBePositive, mustBeFinite}% velocità del suono dell'aria
        temp      {mustBePositive, mustBeFinite}% temperatura dell'aria
        n_vel     = 50;
        n_analisi = 0;
        % Mass
        M         {mustBePositive, mustBeFinite}
        W         {mustBePositive, mustBeFinite}       % [N]
        % Geometry
        D_MR      {mustBePositive, mustBeFinite}       % [m]
        A_D_MR    {mustBePositive, mustBeFinite}       % [m^2]
        D_TR      {mustBePositive, mustBeFinite}       % [m]
        A_D_TR    {mustBePositive, mustBeFinite}       % [m^2]
        b         {mustBePositive, mustBeFinite}       % [m]
        % Aerodynamics
        % tip speeds
        omegaR_MR {mustBePositive, mustBeFinite}       % [m/s]
        omegaR_TR {mustBePositive, mustBeFinite}       % [m/s]
        V_inf     {mustBePositive, mustBeFinite}       % [m/s]
        % Propulsion
        available_power {mustBePositive, mustBeFinite} % [W]
             
    end
    methods
        % compute some mass and geometric property
        function obj = derived_properties(obj)
            % Mass
            obj.M         = obj.fuel_load*obj.M_fuel + obj.M_nofuel;
            obj.W         = obj.M*9.81;                             % [N]
            % Geometry
            obj.D_MR      = obj.R_MR*2;                             % [m]
            obj.A_D_MR    = obj.c_MR*obj.R_MR*obj.N_MR;             % [m^2]
            obj.D_TR      = obj.R_MR*2;                             % [m]
            obj.A_D_TR    = obj.c_TR*obj.R_TR*obj.N_TR;             % [m^2]
            % Distance from the cg to the center of the tail rotor,
            % approximated formula, to be used when b is unknown
            obj.b         = obj.R_MR + obj.R_TR + 0.5;              % [m]
            % Aerodynamics
            % tip speeds
            obj.omegaR_MR = obj.omega_MR;                           % [m/s]
            obj.omegaR_TR = obj.omega_TR;                           % [m/s]
            obj.V_inf     = linspace(0,obj.V_max,obj.n_vel);        % [m/s]
            % Propulsion
            obj.available_power = ...
                obj.engine_number*obj.engine_power;                 % [W]
        end

        % compute ambient conditions
        function obj = ambient(obj)
            [obj.temp, obj.sound_vel, obj.press, obj.rho] = atmosisa(obj.h);
        end

        % compute angular velocity Main Rotor
        function obj = rot_vel_MR(obj,valIN,val)
            switch valIN
                case 'RPM'
                    obj.RPM_MR   = val;
                    obj.n_MR     = obj.RPM_MR/60;
                    obj.omega_MR = obj.n_MR*2*pi;
                case 'n'
                    obj.n_MR    = val;
                    obj.RPM_MR   = obj.n_MR*60;
                    obj.omega_MR = obj.n_MR*2*pi;
                case 'omega'
                    obj.omega_MR = val;
                    obj.n_MR     = obj.omega_MR/(2*pi);
                    obj.RPM_MR   = obj.n_MR*60;
                otherwise
                    mustBeMember(valIN,{'RPM','n','omega'})
            end
        end
            % compute angular velocity Tail Rotor
        function obj = rot_vel_TR(obj,valIN,val)
            switch valIN
                case 'RPM'
                    obj.RPM_TR   = val;
                    obj.n_TR     = obj.RPM_TR/60;
                    obj.omega_TR = obj.n_TR*2*pi;
                case 'n'
                    obj.n_TR    = val;
                    obj.RPM_TR   = obj.n_TR*60;
                    obj.omega_TR = obj.n_TR*2*pi;
                case 'omega'
                    obj.omega_TR = val;
                    obj.n_TR     = obj.omega_TR/(2*pi);
                    obj.RPM_TR   = obj.n_TR*60;
                otherwise
                    mustBeMember(valIN,{'RPM','n','omega'})
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
            % Input:
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