classdef BEMTset_rotor
    properties
        %% Salita assiale
%       flag for Prandtl correction, If 'on', the paramenter B is unused
        P_correction  {mustBeMember(P_correction,{'on','off',})} = 'off';
        %integrals upper bound for the tip correction 
        B       (1,1) {mustBeInRange(B,0,1,'exclude-lower')}     = 1;
        %% Articulated rotor in forward flight
        % Fattore di scorrimento per il calcolo della potenza parassita
        % Pc0 = sigma*Cd_mean/8*(1+k*mu^2)
        k_mu                                                     =3;
        % stop iteration when res < toll
        toll                                                     = 1e-6;
        % Psi domain for the alpha_e map
        Psi     (:,1)                                   = linspace(0,2*pi);
        % fattore di correzione dovuto alla non-uniformita' dell'induzione
        % sul rotore reale: Pc_i = k*lam_i*T_c
        k_i                                                      = 1;
    end
end