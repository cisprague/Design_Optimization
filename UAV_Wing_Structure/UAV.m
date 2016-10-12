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
      NN = length(R);
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
      Iyy = (pi./4.0)*(R(:,2).^4 - R(:,1).^4);
    end
    function F = Sectional_Force(obj, R)
      % Returns the wing's force distribution at a G-force [N].
      % Gravitational acceleration [m/s^2]
      g = obj.Gravity*obj.Gmax;
      % Approximately linear force distribution
      F = g*(obj.Mass/obj.Length)*(1-(obj.Span_Mesh(R)./obj.Length));
    end
    function SD = Spar_Displacement(obj, R)
      % Returns 2 DoF spar displacement at each node. [m, rad].
      SD = CalcBeamDisplacement(obj.Length, obj.E, obj.Second_Moment(R), ...
                                obj.Sectional_Force(R), obj.N_Elements(R));
    end
    function SD = Vertical_Spar_Displacement(obj, R)
      % Returns vertical spar displacement at each node [m].
      SD = obj.Spar_Displacement(R);
      SD = SD(1:2:2*(obj.N_Nodes(R)));
    end
    function Plot_Vertical_Discplacement(obj, R)
      % Plots the spar's spanwise vertical displacement.
      figure;
      plot(obj.Span_Mesh(R), obj.Vertical_Spar_Displacement(R), 'ks-');
      axis equal;
      xlabel('Spanwise Distance [m]');
      ylabel('Vertical Spar Displacement [m]');
    end
    function SS = Spar_Stress(obj, R)
      % Returns the spar's stress at each node [Pa].
      SS = CalcBeamStress(obj.Length, obj.E, R(:,2), ...
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
        r1    = R(I,1);   % First inner radius [m].
        r2    = R(I+1,1); % Second inner radius [m].
        R1    = R(I,2);   % First outer radius [m].
        R2    = R(I+1,2); % Second outer radius [m].
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
      m  = V.*obj.Density;
      ms = VS.*obj.Density;
    end
  end
end
