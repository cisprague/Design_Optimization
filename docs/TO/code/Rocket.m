classdef Rocket
    % This class describes a rocket, modeled as a
    % three dimensional solid in Euclidian space.
    % The rocket has a length (long side), and a
    % width (or diameter), from which its inertia
    % tensor is computed, in order to
    properties
        % State = [Position,Velocity,Quaternion,Euler Rates,Mass]
        % State = [x,y,z,vx,vy,vz,qr,qi,qj,qk,wx,wy,wz,m]
        State;    % Initial state
        Target;   % Target state
        Length;   % Rocket length
        Width;    % Rocket width
        Isp;      % Specific impulse
        NNodes;   % Number of discritization nodes
        StateDim; % Dimension of state
        ContDim;  % Dimension of control
        NAux;     % Number of auxilary variables
        NMain;    % Number of main variables
        NLPDim;   % Number of NPL variables
        LB;       % Lower NLP bounds
        UB;       % Upper NLP bounds
    end

    methods

        function self = Rocket(state,target,L,W,Isp,NNodes)
            self.State    = state;
            self.Target   = target;
            self.Length   = L;
            self.Width    = W;
            self.Isp      = Isp;
            self.NNodes   = NNodes;
            self.StateDim = 14;
            self.ContDim  = 3;
            self.NAux     = 1; % Final time
            self.NMain    = self.StateDim+self.ContDim;
            self.NLPDim   = self.NMain*self.NNodes+self.NAux;
            self.LB       = [];
            self.UB       = [];
            for n=1:self.NNodes
               % Bounds for state per node
               lbs = [-1000;-1000;0;...
                      -100;-100;-100;...
                      0;0;0;0;...
                      -100;-100;-100;...
                      0];
               ubs = [1000;1000;1000;...
                      100;100;100;...
                      1;1;1;1;...
                      100;100;100;...
                      self.State(14)];
               % Bounds for control per node
               lbc     = [0;0;0];
               ubc     = [44000;0;2*pi];
               % Concatenate
               lb      = [lbs;lbc];
               ub      = [ubs;ubc];
               % Concatenate to master
               self.LB = [self.LB;lb];
               self.UB = [self.UB;ub];
            end
            % Final time bounds (auxilary)
            self.LB = [self.LB;10];
            self.UB = [self.UB;7200]; % 2 hours to land
        end

        function PlotTrajectory(self,x)
            [state,control,tf] = self.DecodeNLP(x);
            t = linspace(0,tf,self.NNodes);

            x = state(:,1);
            y = state(:,2);
            z = state(:,3);
            subplot(3,3,1);
            plot3(x,y,z);
            title('Trajectory [m]','Interpreter','latex');
            xlabel('x','Interpreter','latex');
            ylabel('y','Interpreter','latex');
            zlabel('z','Interpreter','latex')
            grid on

            subplot(3,3,2);
            plot(t,x,t,y,t,z);
            title('Position [m]','Interpreter','latex');
            legend('x','y','z');
            legend('Location','best');
            legend('boxoff');

            vx = state(:,4);
            vy = state(:,5);
            vz = state(:,6);
            subplot(3,3,3);
            plot(t,vx,t,vy,t,vz);
            title('Velocity [m/s]','Interpreter','latex');
            legend('v_x','v_y','v_z');
            legend('Location','best');
            legend('boxoff');

            qr = state(:,7);
            qi = state(:,8);
            qj = state(:,9);
            qk = state(:,10);
            subplot(3,3,4);
            plot(t,qr,t,qi,t,qj,t,qk);
            title('Quaternions [ND]','Interpreter','latex');
            legend('q_r','q_i','q_j','q_k');
            legend('Location','best');
            legend('boxoff');

            wx = state(:,11);
            wy = state(:,12);
            wz = state(:,13);
            subplot(3,3,5);
            plot(t,wx,t,wy,t,wz);
            title('Angular Velocity [rad/s]','Interpreter','latex');
            legend('\omega_x','\omega_y','\omega_z');
            legend('Location','best');
            legend('boxoff');
            
            m  = state(:,14);
            subplot(3,3,6);
            plot(t,m);
            title('Mass [kg]','Interpreter','latex');

            thrust = control(:,1);
            subplot(3,3,7);
            plot(t,thrust);
            title('Thrust [N]','Interpreter','latex');

            incl = control(:,2);
            subplot(3,3,8);
            plot(t,incl);
            title('Thrust Inclination [rad]','Interpreter','latex');

            azim = control(:,3);
            subplot(3,3,9);
            plot(t,azim);
            title('Thrust Azimuth [rad]','Interpreter','latex');
        end


        function [x,fval] = Optimise(self,x0,type)
            func    = @self.Objective;
            nvar    = self.NLPDim;
            A       = [];
            b       = [];
            Aeq     = [];
            beq     = [];
            lb      = self.LB;
            ub      = self.UB;
            nonlcon = @self.Constraints;
            if strcmp(type,'fmincon')
                options = optimoptions(@fmincon,'Display','iter-detailed',...
                    'PlotFcn',@optimplotfval,'UseParallel',true,...
                    'ScaleProblem','obj-and-constr',...
                    'Algorithm','sqp',...
                    'MaxFunctionEvaluations',100000);
                [x,fval] = fmincon(func,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
            elseif strcmp(type,'genetic')
                options = optimoptions(@ga,'PlotFcn',@gaplotbestf,...
                    'UseParallel',true,'Display','iter',...
                    'CreationFcn',@gacreationunlinearfeasible);
                [x,fval] = ga(func,self.NLPDim,A,b,...
                              Aeq,beq,lb,ub,nonlcon,options);
            end
        end

        function mass = Objective(self,x)
            % Decode the descision vector
            state = self.DecodeNLP(x);
            % Compute the final mass
            mass  = -state(self.NNodes,14);
        end

        function [c,ceq] = Constraints(self,x)
            % First we decode the descision vector
            [state,control,tf] = self.DecodeNLP(x);
            % We define the time grid
            t   = linspace(0,tf,self.NNodes);
            ceq = [];
            % Trapezoidal quadrature to compute defects
            for n=1:self.NNodes-1;
                s1  = state(n,:).';
                s2  = state(n+1,:).';
                u1  = control(n,:).';
                u2  = control(n+1,:).';
                f1  = self.RHS(s1,u1);
                f2  = self.RHS(s2,u2);
                h1  = t(n+1)-t(n);
                ceq = [ceq;s2-s1-0.5*h1.*(f1+f2)];
            end
            % Satisfy quaternion norms
            for n=1:self.NNodes
                q   = state(n,7:10).';
                ceq =[ceq;q(1)^2+q(2)^+q(3)^2+q(4)^2-1];
            end
            % Satisfy initial condition
            ceq = [ceq;self.State-state(1,:).'];
            % Satisfy terminal conditions (no mass)
            ceq = [ceq;self.Target(1:13)-state(self.NNodes,1:13).'];
            % No inequality constraints
            c   = [];
        end


        function [s,c,tf] = DecodeNLP(self,x)
            tf = x(self.NLPDim); % Final time (auxilary)
            x  = x(1:self.NLPDim-1); % Without time
            x  = reshape(x,[self.NMain,self.NNodes]).';
            s  = x(:,1:self.StateDim);
            c  = x(:,self.StateDim+1:self.NMain);
        end


        function dS = RHS(self,state,control)
            p    = state(1:3);   % Position [m]
            v    = state(4:6);   % Velocity [m]
            q    = state(7:10);  % Quaternion [ND]
            w    = state(11:13); % Angular velocity [rad/s]
            m    = state(14);    % Mass [kg]
            T    = control(1);   % Thrust magnitude [N]
            incl = control(2);   % Thrust body inclination [rad]
            azim = control(3);   % Thrust body azimuth [rad]
            % Thrust direction in body frame
            gamma = [sin(incl)*cos(azim);sin(incl)*sin(azim);cos(incl)];
            ub    = T.*gamma; % Thrust vector in body frame [N]
            u     = self.QuatVecTrans(q,ub); % Inertial thrust vector [N]
            r     = [0;0;-self.Length/2]; % Moment arm
            tau   = cross(r,ub); % Thrust moment [Nm]
            % Inertia of solid
            I     = self.Inertia(self.Width,self.Length,m);
            g     = [0;0;-9.807]; % Gravitational [m/s^2]
            % Assemble the state transition
            dS(1:3,1) = v; %[m/s]
            dS(4:6)   = g + u./m; % Translation [m/s^2]
            dS(7:10)  = self.dqdt(q,w); % Quaternions [ND]
            dS(11:13) = inv(I)*tau; % Rotational [rad/s^2]
            dS(14)    = -T/(self.Isp*9.807); % Mass flow [kg/s]
        end

        function vec = QuatVecTrans(self,q,u)
            qcon  = [q(1);-q(2);-q(3);-q(4)]; % Conjugate
            qnorm = norm(q); % Quaternion normal
            qinv  = qcon./qnorm; % Quaternion inverse
            u     = [0;u]; % Make quaternion
            u     = self.QuatMult(u,qinv);
            u     = self.QuatMult(q,u);
            vec   = u(2:4);
        end

        function qp = QuatMult(self,q1,q2)
            s1 = q1(1);
            v1 = q1(2:4);
            s2 = q2(1);
            v2 = q2(2:4);
            s  = s1*s2 - dot(v1,v2);
            v  = s1.*v2 + s2.*v1 + cross(v1,v2);
            qp = [s;v];
        end

        function I = Inertia(self,r,h,m)
            Ixy = (m/12)*(3*r^2 + h^2);
            Iz = (m*r^2)/2;
            I  = [Ixy,0,0;0,Ixy,0;0,0,Iz];
        end

        function dq = dqdt(self,q,w)
            w  = [0;w]; % Make quaternion
            dq = self.QuatMult(q./2,w);
        end

    end
end
