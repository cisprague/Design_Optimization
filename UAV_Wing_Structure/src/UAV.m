classdef UAV
    properties
        Length;   % Wing semi-span [metres].
        Density;  % Density of spar [kg/m^3].
        E;        % Longitudinal modulus [Pascals].
        Strength; % Ultimate strength [Pa].
        Mass;     % Aircraft operational weight [Kg].
        Gravity;  % Gravity [m/s^2].
        Tmin;     % Minimum spar thickness [m]
        Rimin;    % Minimum inner spar radius [m]
        Romax;    % Maximum outer spar radius [m]
        Gmax;     % Maximum operational G-force.
        NCol;     % Number of collocation points.
        NPert;    % Number of perturbation pts.
    end
    
    methods
        function obj = UAV(L, dens, E, stren, m, ...
                g, tmin, rimin, romax, Gmax,NCol,NPert)
            % Constructs an instance of the UAV class.
            obj.Length   = L;
            obj.Density  = dens;
            obj.E        = E;
            obj.Strength = stren;
            obj.Mass     = m;
            obj.Gravity  = g;
            obj.Tmin     = tmin;
            obj.Rimin    = rimin;
            obj.Romax    = romax;
            obj.Gmax     = Gmax;
            obj.NCol     = NCol;
            obj.NPert    = NPert;
        end
        
        function NN = N_Nodes(obj, R)
            % Returns the integer number of nodes
            % according to the number of radii given.
            NN = length(R)/2;
        end
        
        function NE = N_Elements(obj, R)
            % Returns the integer number of finite
            % elements.
            NE = obj.N_Nodes(R) - 1;
        end
        
        function X = Span_Mesh(obj, R)
            % Returns the spanwise coordinates [m]
            X = linspace(0, obj.Length, obj.N_Nodes(R)).';
        end
        
        function Iyy = Second_Moment(obj, R)
            % Returns the second moment of area for a
            % circular annulus at each spanwise
            % location [m^4].
            for I = 1:obj.N_Nodes(R)
                Iyy(I,1) = (pi/4)*(R(2*I)^4 - R(2*I-1)^4);
            end
        end
        
        function [Fn,Fp,df] = Sectional_Force(obj, R)
            % Returns the wing's force distribution at a
            % G-force [N].
            % Gravitational acceleration [m/s^2]
            g = obj.Gravity*obj.Gmax;
            % Spanwise mesh
            x = obj.Span_Mesh(R);
            % Length of spar
            L = obj.Length;
            % Approximately linear force distribution
            % Nominal force distribution
            Fn = g*(obj.Mass/L)*...
                (1-(x./obj.Length));
            % Initialize probabalistic perturbation
            df = zeros(size(Fn));
            for n=1:4
                % Uniform random variable
                xi = normrnd(0,Fn(1)/(10*n));
                % Probabilistic perturbation
                df = df + xi*cos(((2*n-1)*pi.*x)/(2*L));
            end
            Fp = Fn + df;
        end
        
        function Plot_Sectional_Force(obj, R)
            % Plots the spanwise load distribution.
            [Fn,Fp] = obj.Sectional_Force(R);
            x = obj.Span_Mesh(R)
            figure;
            plot(x,Fn,'ks-',x,Fp,'ko-')
            xlabel('Spanwise Distance [m]');
            ylabel('Sectional Force [N]');
            legend('Nominal','Perturbed')
        end
        
        function [SD, SDV, SDp, SDVp] = Spar_Displacement(obj, R)
            % Nominal and perturbed sectional force
            [F,Fp] = obj.Sectional_Force(R);
            % Returns 2 DoF spar displacement at each
            % node. [m, rad].
            SD = CalcBeamDisplacement(...
                obj.Length, obj.E, obj.Second_Moment(R), ...
                F, obj.N_Elements(R) ...
                );
            SDp = CalcBeamDisplacement(...
                obj.Length, obj.E, obj.Second_Moment(R), ...
                Fp, obj.N_Elements(R) ...
                );
            SDV = SD(1:2:2*(obj.N_Nodes(R)));
            SDVp  = SDp(1:2:2*(obj.N_Nodes(R)));
        end
        
        function Plot_Spar_Displacement(obj, R)
            % Plots the spar's spanwise vertical
            % displacement.
            [SD,SDV,SDp,SDVp] = obj.Spar_Displacement(R);
            x = obj.Span_Mesh(R);
            figure;
            plot(x,SDV,'ks-',x,SDVp,'ko-');
            %axis equal;
            xlabel('Spanwise Distance [m]');
            ylabel('Vertical Spar Displacement [m]');
            legend('Nominal','Perturbed');
        end
        
        function [xi, w] = GaussHermite(obj,n)
            % Computes absciasas x and weights w
            % for Gauss-Hermite quadrature of order n
            i   = 1:n-1;
            a   = sqrt(i/2);
            % Use a diagonal matrix to guarentee real roots
            CM  = diag(a,1) + diag(a,-1);
            % V = column eigenvectors
            % D = diagnonal matrix of eigenvalues
            [V D]   = eig(CM);
            % Get abscissas
            [xi i] = sort(diag(D));
            V       = V(:,i)';
            % Weights
            w       = sqrt(pi) * V(:,1).^2;
        end
        
        function zmax = Max_Height(obj, R)
            % Returns the maximum height at each node [m].
            for I=1:obj.N_Nodes(R)
                zmax(I,1) = R(2*I);
            end
        end
        
        function [SS,SSp] = Spar_Stress(obj, R)
            % Returns the spar's stress at each node [Pa].
            hmax = obj.Max_Height(R);
            [SD,SDV,SDp,SDVp] = obj.Spar_Displacement(R);
            NE = obj.N_Elements(R);
            L = obj.Length;
            E = obj.E;
            SS = CalcBeamStress(L,E,hmax,SD,NE);
            SSp = CalcBeamStress(L,E,hmax,SDp,NE);
        end
        
        function [avg,stdv] = Stochastic_Stress(obj,R)
            % Generate Guass-Hermite quadrature pts.
            [xi,w] = obj.GaussHermite(obj.NCol);
            % We adjust the weights
            w = w./sqrt(pi);
            % Compute the nominal force distribution
            Fn = obj.Sectional_Force(R);
            % Span mesh
            x = obj.Span_Mesh(R);
            % Spar length
            L = obj.Length;
            % Nominal root force
            Fn0 = Fn(1);
            % Second moment
            Iyy = obj.Second_Moment(R);
            % Spar Length
            L = obj.Length;
            % Modulus
            E = obj.E;
            % Number of elements
            NE = obj.N_Elements(R);
            % Number of nodes
            NN = obj.N_Nodes(R);
            % Max Height
            hmax = obj.Max_Height(R);
            % Normal distribution parameters
            mu = 0;
            xir = length(xi);
            sigma = @(n) Fn0/(10*n);
            % Initialize average
            avg = zeros(size(Fn));
            stdv = zeros(size(Fn));
            for i1=1:xir
                pt1 = sqrt(2)*sigma(1)*xi(i1) + mu;
                for i2=1:xir
                    pt2 = sqrt(2)*sigma(2)*xi(i2) + mu;
                    for i3=1:xir
                        pt3 = sqrt(2)*sigma(3)*xi(i3) + mu;
                        for i4=1:xir
                            pt4 = sqrt(2)*sigma(4)*xi(i4) + mu;
                            pts = [pt1;pt2;pt3;pt4];
                            % Probabalistic perturbation
                            df = zeros(size(x));
                            for pti=1:length(pts)
                                df = df + pts(pti).*cos(((2*pti-1)*pi.*x)./(L*2));
                            end
                            % Perturbed force distribution
                            Fp = Fn + df;
                            SD = CalcBeamDisplacement(...
                                obj.Length, obj.E, Iyy, ...
                                Fp, NE);
                            SS = CalcBeamStress(L,E,hmax,SD,NE);
                            weight = w(i1)*w(i2)*w(i3)*w(i4);
                            avg = avg + weight.*SS;
                            stdv = stdv + weight.*SS.^2;
                        end
                    end
                end
            end
            stdv = sqrt(stdv-avg.^2);
        end
        
        function Plot_Stress_Stats(obj,R)
            [mu,stdv] = obj.Stochastic_Stress(R);
            x = obj.Span_Mesh(R);
            figure;
            plot(x,mu,'ks-',x,mu+6.*stdv,'ko-',x,mu-6.*stdv,'ko-');
            xlabel('Spanwise Distance [m]');
            ylabel('Magnitude of Normal Stress [m]');
            legend({'$E(\sigma(x,\xi))$','$E(\sigma(x,\xi)) \pm 6\sqrt{Var(\sigma (x,\xi))}$' },'Interpreter','latex');
        end
        
        function Plot_Spar_Stress(obj, R)
            % Plots the spar's spanwise stress.
            x = obj.Span_Mesh(R);
            [SS,SSp] = obj.Spar_Stress(R);
            figure;
            plot(x,SS,'ks-',x,SSp,'ko-');
            xlabel('Spanwise Distance [m]');
            ylabel('Magnitude of Normal Stress [Pa]');
            legend('Nominal','Perturbed');
        end
        
        function [V, VS] = Spar_Volume(obj, R)
            % Returns volume of spar [m^3].
            % V = Total volume [m^3].
            % VS = Spanwise volume [m^3].
            x = obj.Span_Mesh(R);
            % Loop over all the elements.
            for I=1:obj.N_Elements(R)
                % Index the radii nodes.
                r1    = R(2*I-1,1); % First inner radius [m].
                r2    = R(2*I+1,1); % Second inner radius [m].
                R1    = R(2*I,1);   % First outer radius [m].
                R2    = R(2*I+2,1); % Second outer radius [m].
                x1    = x(I);       % First node location [m].
                x2    = x(I+1);     % Second node location [m].
                % Volume for this element [m^3].
                VS(I,1) = (pi./3)*...
                    (r1.^2-R1.^2+r1*r2+r2.^2-R1*R2-R2.^2)*...
                    (x1-x2);
            end
            V = sum(VS); % Summation of volumes.
        end
        
        function [m, ms] = Spar_Mass(obj, R)
            % Returns the mass of the spar [Kg].
            % m = Total mass [Kg].
            % ms = Spanwise mass [Kg].
            [V, VS] = obj.Spar_Volume(R);
            m       = V.*obj.Density;
            ms      = VS.*obj.Density;
        end
        
        function [lb, ub] = Bound_Constraints(obj, R)
            % Returns the lower and upper bounds for
            % the design vector.
            for I=1:obj.N_Nodes(R)
                lb(2*I-1,1) = obj.Rimin;
                lb(2*I,1)   = obj.Rimin + obj.Tmin;
                ub(2*I-1,1) = obj.Romax - obj.Tmin;
                ub(2*I,1)   = obj.Romax;
            end
        end
        
        function c = Stress_Constraints(obj, R)
            % Returns the nonlinear inequality constraint.
            c = obj.Spar_Stress(R) - obj.Strength;
        end
        
        function [A, b] = Linear_Constraints(obj, R)
            % Returns the linear constraints, dictating
            % that the minimum thickness is abided to.
            for I=1:obj.N_Nodes(R)
                A(I,2*I-1) = 1;
                A(I,2*I)   = -1;
                b(I,1)     = -obj.Tmin;
            end
        end
        
        function jac = Complex_Jacobian(obj, funct, R)
            % Returns the Jacobian matrix of a function
            % using the complex step method.
            fval = funct(R);    % Function value.
            n    = numel(R);    % Number of variables.
            m    = numel(fval); % Number of dependents.
            jac  = zeros(m,n);  % Memory allocation.
            h    = eps(R)*n;    % Small step size.
            for I=1:n
                r        = R;                     % Copy
                r(I,1)   = r(I,1) + h(I,1)*i;     % Perturb
                jac(:,I) = imag(funct(r))/h(I,1); % Diff
            end
            % Transpose to fit dimensions of problem.
            jac = jac.';
        end
        
        function [f, fgrad] = Objective_Function(obj, R)
            f     = obj.Spar_Mass(R);        % Value
            fgrad = obj.Complex_Jacobian(... % Gradient
                @obj.Spar_Mass, R);
        end
        
        function [c, ceq, gradc, gradceq] = Ineq_Constraints(obj, R)
            % Nonlinear constraint
            c       = obj.Stress_Constraints(R);
            ceq     = [];                        % No ineq constraint
            % Gradient of nonlinear constraint
            gradc   = obj.Complex_Jacobian(@obj.Stress_Constraints, R);
            gradceq = []; % No ineq constraint, no gradient.
        end
        
        function c = Stochastic_Stress_Constraints(obj, R)
            % Returns the nonlinear inequality constraint.
            [avg,stdv] = obj.Stochastic_Stress(R);
            c = avg + 6.*stdv - obj.Strength;
        end
        
        function [c,ceq,gradc,gradceq] = Stochastic_Ineq_Constraints(obj,R)
            c = obj.Stochastic_Stress_Constraints(R);
            ceq = [];
            % Complex Jacobian
            gradc = obj.Complex_Jacobian(@obj.Stochastic_Stress_Constraints,R);
            gradceq = [];
        end
        
        function [x,fval,exitflag, ...
                output,lambda,grad,hessian] = Optimize(obj, R, step_type)
            fun      = @obj.Objective_Function; % Objective function.
            x0       = R; % Initial guess of spar radii [m].
            [A, b]   = obj.Linear_Constraints(R);
            Aeq      = [];
            beq      = [];
            [lb, ub] = obj.Bound_Constraints(R);
            nonlcon  = @obj.Stochastic_Ineq_Constraints;
            % Step type specified?
            if ~exist('step_type', 'var')
                step_type = 'central';
            end
            if (step_type == 'central') | (step_type == 'forward')
                % Don't need to specify gradients.
                objgrad = false; constrgrad = false;
            end
            if step_type == 'complex'
                % Must specify gradients.
                objgrad = true; constrgrad = true;
            end
            options  = optimoptions('fmincon', ...
                'Display', 'iter-detailed', ...
                'Algorithm', 'sqp',  ...
                'UseParallel', true, ...
                'PlotFcn', @optimplotfval, ...
                'ScaleProblem', true, ...
                'SpecifyObjectiveGradient', objgrad, ...
                'SpecifyConstraintGradient', constrgrad, ...
                'OptimalityTolerance', 1e-12, ...
                'StepTolerance', 1e-12, ...
                'ConstraintTolerance', 1e-12);
            % Fmincon implementation.
            [x,fval,exitflag, ...
                output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,...
                lb,ub,nonlcon,options)
        end
            
        function Plot_Spar_Shape(obj, R)
            % Returns a plot of the spar's geometry.
            for I=1:obj.N_Nodes(R);
                Ri(I) = R(2*I-1,1);
                Ro(I) = R(2*I,1);
            end
            figure;
            x = obj.Span_Mesh(R);
            plot(x,Ri,'ks-',x,Ro,'ks-',x,-Ri,'ks-',x,-Ro,'ks-');
            xlabel('Spanwise Distance [m]');
            ylabel('Transverse Distance [m]');
        end

    end
end

