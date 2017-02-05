classdef Heat_Exchanger
  properties
    Lx;           % Horizontal width (float)(metres).
    Lymin; Lymax; % Minimum and maximum vertical thicknesses (float)(metres).
    k;            % Thermal conductivity (float)(Watts/metres/Kelvin).
    T1; T2;       % Top and bottom environmental temperatures (float)(Kelvin).
    Nx; Ny;       % Number of elements along x and y directions (integer).
  end
  methods
    function obj = Heat_Exchanger(Lx, Lymin, Lymax, k, T1, T2, Nx, Ny)
      % Constructs Heat Exchanger class instance.
      obj.Lx    = Lx;
      obj.Lymin = Lymin;
      obj.Lymax = Lymax;
      obj.k     = k;
      obj.T1    = T1;
      obj.T2    = T2;
      obj.Nx    = Nx;
      obj.Ny    = Ny;
    end
    function ts = Surface(obj, a)
      % Generates the form of the top surface as a function of (a,x),
      % according to the specified number of coefficients.
      % First necessary coefficient.
      eqn   = @(x) a(1);
      % Summation of the sin series according to length of vector.
      for n = 2:length(a);
        eqn = @(x) eqn(x) + a(n)*sin((2*pi*(n-1)*x)/obj.Lx);
      end
      ts    = eqn(linspace(0, obj.Lx, obj.Nx + 1)).';
    end
    function f = Neg_Flux(obj,a)
      % Calculates the negative flux of the given geometry.
      % Generate top surface mesh.
      h    = obj.Surface(a);
      % Calculate the flux per unit length.
      flux = CalcFlux(obj.Lx, h, obj.Nx, obj.Ny, obj.k, obj.T2, obj.T1);
      % Negative sign to convert for minimization.
      f    = -flux;
    end
    function [c, ceq] = Thickness_Limit(obj, a)
      % Mesh of top surface.
      ts   = obj.Surface(a);
      % Maximum point.
      tmax = max(ts);
      % Minimum point.
      tmin = min(ts);
      % Constraint vector.
      c    = [tmax - obj.Lymax; -tmin + obj.Lymin];
      ceq  = [];
    end
    function [aopt, fval] = Optimise(obj, a0, alg)
      % Specify the objective function to be minimized.
      fun     = @obj.Neg_Flux;
      % Specify the nonlinear inequality constraints.
      nonlcon = @obj.Thickness_Limit;
      % Ignore other specifications.
      A   = []; b   = []; Aeq = []; beq = []; lb  = []; ub  = [];
      % Specify optimization options.
      options = optimoptions(           ...
        'fmincon',                      ... % The optimisation algorithm.
        'Display',              'iter', ... % Display the optimisation output.
        'UseParallel',            true, ... % Use multiple CPUs.
        'Algorithm',               alg, ... % The specified algorithm.
        'MaxFunctionEvaluations', 5000, ... % Many evaluations...
        'MaxIterations',          5000  ... % Many iterations...
      );
      [aopt, fval] = fmincon(fun,a0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    end
    function Visualise(obj, aopt)
      ts   = obj.Surface(aopt);
      bs   = zeros(obj.Nx + 1);
      x    = linspace(0, obj.Lx, obj.Nx + 1);
      area(x, ts);
      tstr = 'Heat Exchanger Geometry ($$Flux = ';
      tstr = strcat(tstr, num2str(-obj.Neg_Flux(aopt)), '~\frac{W}{m}$$)');
      title(tstr, 'Interpreter', 'latex');
      xlabel('Width $$L_x$$ [metres]', 'Interpreter', 'latex')
      ystr = 'Thickness $$L_y$$ [metres]'
      ystr = strcat(ystr, '$$~\mathbf{a} = ', mat2str(aopt), '$$')
      ylabel(ystr, 'Interpreter', 'latex')
    end
  end
end
