classdef BEMTset_rotor
    properties
        % flag for Prandtl correction
        P_correction  {mustBeMember(P_correction,{'on','off',})} = 'off';
        % integrals upper bound for the Pradtl Correction 
        B       (1,1) {mustBeInRange(B,0,1,'exclude-lower')}     = 1
        % stop iteration when res < toll
        toll                                                     = 1e-6;
        % Psi domain for the alpha_e map
        Psi     (:,1)                                   = linspace(0,2*pi);
    end
end