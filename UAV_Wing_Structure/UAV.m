classdef UAV
  properties
    Length;   % Wing semi-span [metres].
    Density;  % Density of spar [kg/m^3].
    E;        % Longitudinal elastic modulus [Pascals].
    Strength; % Ultimate tensile/compresive strength [Pa].
    Mass;     % Aircraft operational weight [Kg].
    Gravity;  % Gravitational acceleration [m/s^2].
    Tmin;     % Minimum thickenss of spar [m]
    Rimin;    % Minimum inner spar radius [m]
    Romax;    % Maximum outer spar radius [m]
    Gmax;     % Maximum operational G-force.
  end

  methods
    function obj = UAV(L, dens, E, stren, m, g, tmin, rimin, romax, Gmax)
      % Constructs an instance of the UAV class.
      obj.Length   = L;     obj.Density = dens;  obj.E       = E;
      obj.Strength = stren; obj.Mass    = m;     obj.Gravity = g;
      obj.Tmin     = tmin;  obj.Rimin   = rimin; obj.Romax   = romax;
      obj.Gmax     = Gmax;
    end

    function NN = N_Nodes(obj, R)
      % Returns the integer number of nodes according
      % to the number of radii given.
      NN = length(R)/2;
    end

    function NE = N_Elements(obj, R)
      % Returns the integer number of finite elements.
      NE = obj.N_Nodes(R) - 1;
    end

    function X = Span_Mesh(obj, R)
      % Returns the spanwise coordinates [m]
      X = linspace(0, obj.Length, obj.N_Nodes(R)).';
    end

    function Iyy = Second_Moment(obj, R)
      % Returns the second moment of area for a circular
      % annulus at each spanwise location [m^4].
      for I = 1:obj.N_Nodes(R)
        Iyy(I,1) = (pi/4)*(R(2*I)^4 - R(2*I-1)^4);
      end
    end

    function F = Sectional_Force(obj, R)
      % Returns the wing's force distribution at a G-force [N].
      % Gravitational acceleration [m/s^2]
      g = obj.Gravity*obj.Gmax;
      % Approximately linear force distribution
      F = g*(obj.Mass/obj.Length)*(1-(obj.Span_Mesh(R)./obj.Length));
    end

    function [SD, SDV] = Spar_Displacement(obj, R)
      % Returns 2 DoF spar displacement at each node. [m, rad].
      SD = CalcBeamDisplacement(obj.Length, obj.E, obj.Second_Moment(R), ...
                                obj.Sectional_Force(R), obj.N_Elements(R));
      SDV = SD(1:2:2*(obj.N_Nodes(R)));
    end

    function Plot_Spar_Displacement(obj, R)
      % Plots the spar's spanwise vertical displacement.
      [SD, SDV] = obj.Spar_Displacement(R);
      figure;
      plot(obj.Span_Mesh(R), SDV, 'ks-');
      axis equal;
      xlabel('Spanwise Distance [m]');
      ylabel('Vertical Spar Displacement [m]');
    end

    function zmax = Max_Height(obj, R)
      for I=1:obj.N_Nodes(R)
          zmax(I,1) = R(2*I);
      end
    end

    function SS = Spar_Stress(obj, R)
      % Returns the spar's stress at each node [Pa].
      SS = CalcBeamStress(obj.Length, obj.E, obj.Max_Height(R), ...
                          obj.Spar_Displacement(R), obj.N_Elements(R));
    end

    function Plot_Spar_Stress(obj, R)
      % Plots the spar's spanwise stress.
      figure;
      plot(obj.Span_Mesh(R), obj.Spar_Stress(R), 'ks-');
      xlabel('Spanwise Distance [m]');
      ylabel('Magnitude of Normal Stress [Pa]')
    end

    function [V, VS] = Spar_Volume(obj, R)
      % Returns volume of spar [m^3].
      % V = Total volume [m^3], VS = Spanwise volume [m^3].
      x = obj.Span_Mesh(R);
      % Loop over all the elements.
      for I=1:obj.N_Elements(R)
        % Index the radii nodes.
        r1    = R(2*I-1,1);   % First inner radius [m].
        r2    = R(2*I+1,1); % Second inner radius [m].
        R1    = R(2*I,1);   % First outer radius [m].
        R2    = R(2*I+2,1); % Second outer radius [m].
        x1    = x(I);     % First node location [m].
        x2    = x(I+1);   % Second node location [m].
        % Volume for this element [m^3].
        VS(I,1) = (pi./3)*(r1.^2-R1.^2+r1*r2+r2.^2-R1*R2-R2.^2)*(x1-x2);
      end

      V = sum(VS);
    end

    function [m, ms] = Spar_Mass(obj, R)
      % Returns the mass of the spar [Kg].
      % m = Total mass [Kg], ms = Spanwise mass [Kg].
      [V, VS] = obj.Spar_Volume(R);
      m       = V.*obj.Density;
      ms      = VS.*obj.Density;
    end

    function [lb, ub] = Bound_Constraints(obj, R)
      for I=1:obj.N_Nodes(R)
        lb(2*I-1,1) = obj.Rimin;
        lb(2*I,1)   = obj.Rimin + obj.Tmin;
        ub(2*I-1,1) = obj.Romax - obj.Tmin;
        ub(2*I,1)   = obj.Romax;
      end
    end

    function c = Stress_Constraints(obj, R)
      c = obj.Spar_Stress(R) - obj.Strength;
    end

    function [A, b] = Linear_Constraints(obj, R)
      for I=1:obj.N_Nodes(R)
        A(I,2*I-1) = 1;
        A(I,2*I)   = -1;
        b(I,1)     = -obj.Tmin;
      end
    end

    function jac = Complex_Jacobian(obj, funct, R)
      fval = funct(R);
      n    = numel(R);
      m    = numel(fval);
      jac  = zeros(m,n);
      h    = eps(R)*n;
      for I=1:n
        r        = R;
        r(I,1)   = r(I,1) + h(I,1)*i;
        jac(:,I) = imag(funct(r))/h(I,1);
      end
      jac = jac.';
    end

    function [f, fgrad] = Objective_Function(obj, R)
      f     = obj.Spar_Mass(R);
      fgrad = obj.Complex_Jacobian(@obj.Spar_Mass, R);
    end

    function [c, ceq, gradc, gradceq] = Ineq_Constraints(obj, R)
      c       = obj.Stress_Constraints(R);
      ceq     = [];
      gradc   = obj.Complex_Jacobian(@obj.Stress_Constraints, R);
      gradceq = [];
    end

    function [x,fval] = Optimize(obj, R, step_type)
      fun      = @obj.Objective_Function; % Objective function.
      x0       = R; % Initial guess of spar radii [m].
      [A, b]   = obj.Linear_Constraints(R);
      Aeq      = [];
      beq      = [];
      [lb, ub] = obj.Bound_Constraints(R);
      nonlcon  = @obj.Ineq_Constraints;
      if ~exist('step_type', 'var')
          step_type = 'central';
      end
      if (step_type == 'central') | (step_type == 'forward')
        options  = optimoptions('fmincon', ...
                        'Display', 'iter-detailed', ...
                        'Algorithm', 'sqp',  ...
                        'UseParallel', true, ...
                        'PlotFcn', @optimplotfval, ...
                        'ScaleProblem', true, ...
                        'OptimalityTolerance', 1e-10, ...
                        'FiniteDifferenceType', step_type);
      end
      if step_type == 'complex'
        options  = optimoptions('fmincon', ...
                        'Display', 'iter-detailed', ...
                        'Algorithm', 'sqp',  ...
                        'UseParallel', true, ...
                        'PlotFcn', @optimplotfval, ...
                        'ScaleProblem', 'obj-and-constr', ...
                        'OptimalityTolerance', 1e-15, ...
                        'SpecifyObjectiveGradient', true, ...
                        'SpecifyConstraintGradient', true);
      end
      % Fmincon implementation.
      [x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    end

    function Plot_Spar_Shape(obj, R)
      for I=1:obj.N_Nodes(R)
        Ri(I) = R(2*I-1,1);
        Ro(I) = R(2*I,1);
      end
      figure;
      x = obj.Span_Mesh(R);
      plot(x,Ri,x,Ro,x,-Ri,x,-Ro);
    end

  end
end
