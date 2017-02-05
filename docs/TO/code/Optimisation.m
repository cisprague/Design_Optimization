clear all
% Instantiate the wing
Raptor = UAV(... % With properties:
  7.5,       ... % Wing semi-span [m].
  1600,      ... % Spar density [kg/m^3].
  7e10,      ... % Modulus of elasticity [Pa].
  600e6,     ... % Maximum spar strength [Pa].
  500,       ... % Mass of aircraft [Kg].
  9.807,     ... % Earth's gravity [m/s^2].
  2.5e-3,    ... % Minimum spar thickness [m].
  1e-2,      ... % Minimum inner spar radius [m].
  5e-2,      ... % Maximum outer spar radius [m].
  2.5,       ... % Maximum operational G-force.
  5,         ... % Number of collocation points.
  4          ... % Number of perturbations pts.
  );

% Take an initial guess of the spar's radii.
R = zeros(80,1);
for I=1:Raptor.N_Nodes(R)
    R(2*I-1,1) = 0.0415;
    R(2*I, 1)  = 0.05;
end
[x,fval,exitflag,output,lambda,grad,hessian] = Raptor.Optimize(R,'complex')
Raptor.Plot_Stress_Stats(x)