%% Heat Exchanger
% Authored by Christopher Iliffe Sprague
% Christopher.Iliffe.Sprague@gmail.com
% https://github.com/CISprague/Design-Optimization.git

%% Heat Exchanger Class
classdef Heat_Exchanger
    % Defines a heat exchanger with straight bottom, left, and right sides,
    % as well as a user designed top surface that is variable in the horizontal
    % direction. The heat exchanger's top and bottom temperatures, thermal
    % conductivity, and optimization resolutions are provided by the user.

    %% Properties
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
        % Note: Ly must be defined as an anonymous function with
        % variables (a,x), where a is a vector of coefficients of
        % user define length, and x is the horizontal distance
        % along the heat exchanger.
    end

    %% Methods
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
            % Generate top surface mesh
            h    = obj.Top_Surface(a);
            % Calculate the flux per unit length
            flux = CalcFlux(obj.Lx, h,obj.Nx,obj.Ny,obj.k,obj.T2,obj.T1);
            % Negative sign to convert maximization to minimization
            f    = -flux;
        end
        function [c, ceq] = thickness_limit(obj, a)
            % Mesh of top surface
            ts   = obj.Top_Surface(a);
            % Heighest point
            tmax = max(ts);
            % Minimum point
            tmin = min(ts);
            % Constraint vector
            c    = [tmax - obj.Lymax; -tmin + obj.Lymin];
            ceq  = [];
        end
        function [aopt, fval] = Optimize(obj, a0)
            fun     = @obj.objfun;
            nonlcon = @obj.thickness_limit;
            A       = [];
            b       = [];
            Aeq     = [];
            beq     = [];
            lb      = [];
            ub      = [];
            options = optimoptions('fmincon','Display', 'iter');
            [aopt, fval] = fmincon(fun,a0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        end
    end
end
