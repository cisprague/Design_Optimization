%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heat Exchanger                                       %
% Authored by Christopher Iliffe Sprague               %
% Christopher.Iliffe.Sprague@gmail.com                 %
% https://github.com/CISprague/Design-Optimization.git %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Heat_Exchanger
    properties
        % Horrizontal width and vertical thickness (cm)
        Lx; Ly
        % Minimum and maximum vertical thicknesses (cm)
        Lymin; Lymax
        % Thermal conductivity (W/mK)
        k
        % Environmental temperatures (c)
        T1; T2
    end
    methods
        function Heat_Flux(obj)
            DT = obj.T1 - obj.T2;
            Int = HE.Ly
        end
    end
end
