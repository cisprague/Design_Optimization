classdef Car
    properties
        % Fixed parameters
        r1; alpha0; Q0;
        % Design parameters
        omega; r2; alpha1;
        % Bounds
        omegab; r2b; alpha1b;
        % Extra
        T; samps; noise;
    end
    methods
        function self = Car(      ...
            r1, alpha0, Q0,       ...
            omega, r2, alpha1,    ...
            omegab, r2b, alpha1b, ...
            T, samps, noise       ...
        )
            % Fixed parameters
            self.r1      = r1;
            self.alpha0  = alpha0;
            self.Q0      = Q0;

            % Design parameters
            self.omega   = omega;
            self.r2      = r2;
            self.alpha1  = alpha1;

            % Bounds
            self.omegab  = omegab;
            self.r2b     = r2b;
            self.alpha1b = alpha1b;

            % Integration period
            self.T       = T;
            % # True function samples
            self.samps   = samps;
            % Amount of noise
            self.noise   = noise;
        end
        function D = Design_Vector(self)
            D = [self.omega, ...
                 self.r2,    ...
                 self.alpha1];
        end
        function B = Bounds(self)
          B = [self.omegab; ...
               self.r2b;    ...
               self.alpha1b];
        end
        function dphi = Rate(self, tau, phi)
            % Described the Tilt-A-Whirl's
            % second order ordinary differential
            % equations of motion.
            p         = phi(1); % p
            dp        = phi(2); % p'
            dphi(1,1) = dp; % p''
            % The system dynamics from
            % Eq. (27) of "Chaos at the amusement
            % park: Dynamics of the Tilt-A-Whirl
            % (Kautz and Huggard):
            g     = 9.802;
            gamma = 1/(3*self.omega);
            gamma = sqrt(g/self.r2)*gamma;
            eps   = self.r1/(9*self.r2);
            alpha = -self.alpha1*cos(tau);
            alpha = alpha+self.alpha0;
            beta  = 3*self.alpha1*sin(tau);
            d2p   = -(eps-(gamma^2)*alpha)*sin(p);
            d2p   = d2p-(gamma/self.Q0)*dp;
            d2p   = d2p-gamma^2*beta*cos(p);
            dphi(2,1) = d2p;
        end
        function [phi, t] = Propogate(self)
            % Propogates dynamics over over
            % dimensionless time.
            T        = self.T;
            [t, phi] = ode45(@self.Rate,[0,T],[0,0]);
        end
        function [stdp, phi, t] = StdPhi(self)
            % Returns standard deviations of dphi/dt.
            [phi,t] = self.Propogate();
            % Sample terminal dynamics only.
            el      = length(phi);  % End index
            tl      = round(el/3);  % Term index
            y       = phi(tl:el,2); % Term sample
            yb      = mean(y);
            stdp    = mean((y-yb).^2);
            stdp    = 3*self.omega*sqrt(stdp);
        end
        function [fval,s2] = GP_Surrogate(self, x)
          [hyp, xs, ys] = self.GP_Hyper();
          % Evaluate mean at variables
          [fval, s2] = gp(                 ...
            hyp, @infExact, [], @covSEiso, ...
            @likGauss, xs, ys, x           ...
          );
        end
        function [x,fval] = Optimize(self, alg)
          fun = @(x) -self.GP_Surrogate(x);
          % Initial guess
          x0  = self.Design_Vector();
          nv = length(x0);
          % Lower bounds
          lb  = [self.omegab(1), ...
                 self.r2b(1),    ...
                 self.alpha1b(1)];
          % Upper bounds
          ub  = [self.omegab(2), ...
                 self.r2b(2),    ...
                 self.alpha1b(2)];
          % Extraneous constraints
          A   = []; b   = [];
          Aeq = []; beq = []; nonlcon = [];
          % Optimize the objective
          disp('Optimising Guassian surrogate...')

          if strcmp(alg,'simulated-annealing')
            options = optimoptions(@simulannealbnd, ...
              'PlotFcn',{@saplotbestf,@saplottemperature, ...
                         @saplotf,@saplotstopping});
            [x,fval] = simulannealbnd(fun,x0,lb,ub, options);

          elseif strcmp(alg,'particle-swarm')
            options = optimoptions(@particleswarm, ...
              'UseParallel', true, ...
              'HybridFcn', @fmincon, ...
              'PlotFcn', @pswplotbestf, ...
              'Display', 'iter' ...
            );
            [x,fval] = particleswarm(fun,nv,lb,ub, options);

          elseif strcmp(alg,'genetic')
            options = optimoptions(@ga, ...
              'PlotFcn', @gaplotbestf, ...
              'Display', 'iter', ...
              'UseParallel', true ...
            );
            [x,fval] = ga(fun,nv,A,b,Aeq,beq,lb,ub,nonlcon,options);
          elseif strcmp(alg,'fmincon')
            options = optimoptions(@fmincon, ...
              'Display', 'iter', ...
              'PlotFcn', {@optimplotfval,@optimplotx}, ...
              'UseParallel', true, ...
              'ScaleProblem', 'obj-and-constr' ...
            );
            [x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

          end
        end

        function [hyp, xs, ys] = GP_Hyper(self)
          % Hyperparameters exist?
          if exist('Hyper.mat','file') == 2
            load('Hyper.mat');
            % Want a different sampling?
            if (self.samps ~= length(xs) ...
                || self.noise ~= ns      ...
                || self.T ~= T)
              disp('Different sampling...');
              need = true; % Then sample again
            else
              need = false; % Otherwise don't
            end
          else
            disp('Need sampling...');
            need = true; % No such file :(
          end

          if need == true
            % Latin hyper cube sampling
            [ys, xs]  = self.LHS_Sampling();
            % Set Matern covariance function
            covfunc   = {@covMaterniso, 1};
            % hyp.cov = [log(1); lfilefileog(sigma)];
            hyp.cov   = log([1/4;1.0]);
            % Gaussian likelihood function
            likfunc   = @likGauss;
            % Noise
            ns = self.noise;
            hyp.lik   = log(ns);
            % Minimise negative log
            % likelihood function
            disp('Minimising negative log likelihood function...');
            hyp = minimize(...
              hyp, @gp, -100, @infExact, ...
              [], @covSEiso, @likGauss, xs, ys ...
            );
            T = self.T;
            % Save for later
            save('Hyper','hyp','xs','ys','ns','T')
          end
        end
        function [stphi, x] = LHS_Sampling(self)
          % Define sample points with latin
          % hyper cube sampling.
          % # Variables
          n = length(self.Design_Vector());
          s = self.samps; % # of samples
          x = lhsdesign(s,n);
          % Variable bounds
          xl = [self.omegab(1);
                self.r2b(1);
                self.alpha1b(1)];
          xu = [self.omegab(2);
                self.r2b(2);
                self.alpha1b(2)];
          % Scaling according to bounds
          for i = 1:n
            x(:,i) = x(:,i).*(xu(i) - xl(i));
            x(:,i) = x(:,i) + xl(i);
          end
          % Evaluate objective function
          % at LHS samples
          ndes = length(x); % # designs
          for i = 1:ndes
            self.omega  = x(i,1);
            self.r2     = x(i,2);
            self.alpha1 = x(i,3);
            stphi(i,1)  = self.StdPhi();
          end
        end
        function [X,Y,F] = Surface(self, cvar, np)
          % The constant variable index
          if strcmp(cvar,'omega')
            ci = 1;
          elseif strcmp(cvar,'r2')
            ci = 2;
          elseif strcmp(cvar,'alpha1')
            ci = 3;
          end
          % Design variable bounds
          B     = self.Bounds();
          ind   = [1,2,3];
          % Indicies of variables
          ind   = ind(ind~=ci);
          ix    = ind(1);
          iy    = ind(2);
          % Bounds of variables
          B     = [B(ix,:); B(iy,:)];
          % Axes labels
          axst  = {'$\omega$', '$r_2$', '$\alpha_1$'};
          axstx = axst(ix);
          axsty = axst(iy);
          % Generate grid
          x     = linspace(B(1,1),B(1,2),np);
          y     = linspace(B(2,1),B(2,2),np);
          [X,Y] = meshgrid(x,y);
          % Constant variable
          vars = self.Design_Vector();
          z = vars(ci);
          % Objective values
          dim = size(X);
          for j = 1:dim(2)
            x = X(:,j);
            y = Y(:,j);
            z = ones(size(x)).*z;
            d(:,ix) = x;
            d(:,iy) = y;
            d(:,ci) = z;
            F(:,j) = self.GP_Surrogate(d);
          end
          surfc(X,Y,F);
          xlabel(axstx,'Interpreter','LaTex');
          ylabel(axsty,'Interpreter','LaTex');
          zlabel('$\sigma(\frac{d\phi}{dt})$', ...
            'Interpreter', 'LaTex')
        end
    end
end
