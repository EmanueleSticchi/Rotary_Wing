classdef BEMTset
    properties
        toll   (1,1){mustBePositive,mustBeFinite}=1e-6; % Tollerance
        P_correction  {mustBeMember(P_correction,{'on','off',})} = 'off';
        Hub_correction {mustBeMember(Hub_correction,{'on','off',})} = 'off';
        Cd_hub (1,1) {mustBeNumeric,mustBeFinite}=0.5;
        Design {mustBeMember(Design,{'on','off',})} = 'off';
    end
end

