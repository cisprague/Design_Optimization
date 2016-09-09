%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heat Exchanger                                       %
% Authored by Christopher Iliffe Sprague               %
% Christopher.Iliffe.Sprague@gmail.com                 %
% https://github.com/CISprague/Design-Optimization.git %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Heat_Exchanger
    properties
        % Horrizontal width and vertical thickness (cm)
        Lx; Ly;
        % Minimum and maximum vertical thicknesses (cm)
        Lymin; Lymax;
        % Thermal conductivity (W/mK)
        k;
        % Environmental temperatures (c)
        T1; T2;
        % Number of elements along x direction and y direction
        Nx; Ny
        % Horizontal domain
    end
    methods
        function h = Top_Surface(obj, a)
            % Insert coefficients
            h = @(x) obj.Ly(a, x);
            % Horizontal domain
            x = linspace(0, obj.Lx, obj.Nx+1);
            % Top surface mesh
            h = h(x).';
        end
        function f = objfun(obj,a)
            % Calculate flux per unit length
            h = obj.Top_Surface(a);
            flux = CalcFlux(obj.Lx, h,obj.Nx,obj.Ny,obj.k,obj.T2,obj.T1);
            f = -flux;
        end
        function [c, ceq] = confun(obj, a)
            hmax = max(obj.Top_Surface(a));
            hmin = min(obj.Top_Surface(a));
            c = [hmax - obj.Lymax; hmin + obj.Lymin];
            ceq = [];
        end
        function opt = Optimize(obj, a0)
            [aopt,fval] = fmincon(@objfun,a0,[],[],[],[],[],[],@confun);
            opt = aopt
        end
    end
end