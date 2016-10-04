% Instantiate the wing
Raptor = UAV(...
  7.5, ... % Wing semi-span [m]
  7e10, ... % Modulus of elasticity [Pa]
  10, ... % Number of finite elements
  500)

F =Raptor.Sectional_Force()
plot(F)
