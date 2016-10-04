classdef UAV
  properties
    Length; % Wing semi-span [metres]
    E;      % Longitudinal elastic modulus [Pascals]
    N_Elem; % Number of finite elements to use
    Mass;   % Aircraft operational weight [Kg]
  end
  methods
    function obj = UAV(L, E, NE, m)
      % Givens
      obj.Length = L;
      obj.E      = E;
      obj.N_Elem = NE;
      obj.Mass   = m;
    end
    function F = Sectional_Force(obj)
      x = linspace(0, obj.Length, obj.N_Elem);
      F = (obj.Mass/obj.Length)*(1-(x/obj.Length))
    end
  end
end
