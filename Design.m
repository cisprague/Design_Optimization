% Implementation
% Authored by Christopher Iliffe Sprague
% Christopher.Iliffe.Sprague@gmail.com
% https://github.com/CISprague/Design-Optimization.git

clear

%% Givens
% Instantiate the heat exchanger.
HE = Heat_Exchanger;
% Horizontal width and height (cm).
HE.Lx = 5;
% Vertical thickness (cm)
HE.Ly = @(a,x) a(1) + a(2)*sin((2*pi*(2-1)*x)/HE.Lx) + a(3)*sin((2*pi*(3-1)*x)/HE.Lx) + a(4)*sin((2*pi*(4-1)*x)/HE.Lx);
% Minimum and maximum thickness constraints (cm).
HE.Lymin = 1; HE.Lymax = 5;
% Thermal conductivity (w/m/k).
HE.k = 20;
% Temperature of water and air, respectively (c).
HE.T1 = 90; HE.T2 = 20;
% Number of elements along x direction and y direction, resepectively.
HE.Nx = 100; HE.Ny = 100;

%% Implementation
% Initial guess within contraints
a0 = [4;.75;0.2;0];
[aopt, fval] = HE.Optimize(a0)
x = linspace(0, HE.Lx, HE.Nx + 1);
plot(x.', HE.Top_Surface(aopt))