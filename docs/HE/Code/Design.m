%% Heat Exchanger Optimisation
% Authored by Christopher Iliffe Sprague
% Christopher.Iliffe.Sprague@gmail.com
% https://github.com/CISprague/Design-Optimization.git

% Clear the workspace.
clear
%sqp
aopt_sqp = [0.0300; -0.0019; -0.0015; 0.0008; -0.0190];
%
aopt_ipopt = [0.0300; -0.0025; 0.0044; 0.0032; 0.0145; -0.0048];

% Instantiate a Heat Exchanger object with the constructor.
HE = Heat_Exchanger(...
  0.05,             ... % Horizontal width (metres).
  0.01,             ... % Minimum thickness (metres).
  0.05,             ... % Maximim thickness (metres).
  20,               ... % Thermal conductivity (Watts•metres⁻¹•Kelvin⁻¹).
  90 + 273.15,      ... % Temperature of bottom surface (Kelvin).
  20 + 273.15,      ... % Temperature of top surface (Kelvin).
  200,              ... % Number of horizontal elements.
  200               ... % Number of vertical elements.
);

HE.Visualize(aopt_sqp)
