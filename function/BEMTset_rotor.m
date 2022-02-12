classdef BEMTset_rotor
    properties
        P_correction  {mustBeMember(P_correction,{'on','off',})} = 'off';
        Design {mustBeMember(Design,{'on','off',})} = 'off';
    end
end