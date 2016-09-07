%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation                                       %
% Authored by Christopher Iliffe Sprague               %
% Christopher.Iliffe.Sprague@gmail.com                 %
% https://github.com/CISprague/Design-Optimization.git %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%% Instantiate the Heat Exchanger
HE = Heat_Exchanger;
% Set width (cm)
HE.Lx = 5;
% Polynomial thickenss (cm)
HE.Ly = @(a1, a2, a3, a4, x) a1 + a2*x + a3*x.^2 + a4*x.^3;
% Minimum and maximum thicknesses (cm)
HE.Lymin = 1; HE.Lymax = 5;
% Thermal conductivity (W/mK)
HE.k = 20;
% Environmental temperatures (c)
HE.T1 = 90; HE.T2 = 20;

% Example
HE.Ly = @(x) HE.Ly(1,2,3,4,x)