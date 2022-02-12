classdef BEMTset_rotor
    properties
        P_correction  {mustBeMember(P_correction,{'on','off',})} = 'off';
        B       (1,1) {mustBeInRange(B,0,1,'exclude-lower')} = 1
    end
end