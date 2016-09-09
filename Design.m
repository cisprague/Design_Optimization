%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation                                       %
% Authored by Christopher Iliffe Sprague               %
% Christopher.Iliffe.Sprague@gmail.com                 %
% https://github.com/CISprague/Design-Optimization.git %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%% Instantiate heat exchanger
HE = Heat_Exchanger;
% Horizontal width (cm) and vertical thickness as a function of x (cm)
HE.Lx = 5; HE.Ly = @(a, x) a(1) + a(2)*x + a(3)*x.^2;
% Minimum and maximum thickness (cm)
HE.Lymin = 1; HE.Lymax = 5;
% Thermal conductivity (w/m/k)
HE.k = 20;
% Temperature of water and air (c)
HE.T1 = 90; HE.T2 = 20;
% Number of elements along x direction and y direction
HE.Nx = 100; HE.Ny = 100;

%% Optimization
% Test the objective function
a0 = [1;2;3];
HE.objfun(a0)
HE.confun(a0);
%aopt = HE.Optimize(a0)
%plot(HE.Top_Surface(aopt))