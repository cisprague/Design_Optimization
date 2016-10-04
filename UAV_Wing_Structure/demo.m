% demonstrate the CalcBeamDisplacement and CalcBeamStress routines
clear all;
close all;

L = 1.0;
E = 1.0;

Nelem = 10;
Iyy = ones(Nelem+1,1);
force = ones(Nelem+1,1);

[u] = CalcBeamDisplacement(L, E, Iyy, force, Nelem);

% plot the vertical displacement
x = [0:L/Nelem:L];
plot(x,u(1:2:2*(Nelem+1)),'ks-');
xlabel('distance along wing')
ylabel('vertical displacement of spar')

% plot the stresses
figure;
zmax = ones(Nelem+1,1);
[sigma] = CalcBeamStress(L, E, zmax, u, Nelem);
plot(x,sigma,'ks-')
xlabel('distance along wing')
ylabel('magnitude of normal stress')