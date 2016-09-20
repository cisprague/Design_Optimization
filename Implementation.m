%% Heat Exchanger Optimisation
% Authored by Christopher Iliffe Sprague
% Christopher.Iliffe.Sprague@gmail.com
% https://github.com/CISprague/Design_Optimization

% Clear the workspace.
clear

% Instantiate a Heat Exchanger object with the constructor.
HE = Heat_Exchanger(...
  0.05,             ... % Horizontal width (metres).
  0.01,             ... % Minimum thickness (metres).
  0.05,             ... % Maximim thickness (metres).
  20,               ... % Thermal conductivity (Watts•metres⁻¹•Kelvin⁻¹$).
  90 + 273.15,      ... % Temperature of bottom surface (Kelvin).
  20 + 273.15,      ... % Temperature of top surface (Kelvin).
  200,              ... % Number of horizontal elements.
  200               ... % Number of vertical elements.
);

% Initial guess.
a0 = [0.03; 0.0; 0.0; 0.0; 0.0];
% Optimisation.
[aopt, fval] = HE.Optimise(a0);
% Visualization.
HE.Visualise(aopt);

% aopt_sqp_5_100   = 'infeasible'
% aopt_sqp_6_100   = 'infeasible';
% aopt_sqp_5_200   = 'infeasible';
% aopt_sqp_6_200   = [0.0300; -0.0019; -0.0015; 0.0008; -0.0190; 0.0001];
% aopt_ipopt_5_100 = 'infeasible';
% aopt_ipopt_6_100 = [0.0300; -0.0025; 0.0044; 0.0032; 0.0145; -0.0048];
% aopt_ipopt_5_200 = [0.0300; -0.0021; 0.0024; -0.0187; -0.0025];
% aopt_ipopt_6_200 = [0.0300; -0.0025; 0.0044; 0.0032; 0.0146; -0.0049];
