f = @(x) (x-3).*x.^3.*(x-6).^4;

xi = [-2.351; -1.336; -.4361; 2.351; 1.336; .4361];
omega = [.00453; .1571; .7245; -.00453; -.1571; -.7245];

m = 6;
sigma =1;

% Scaling
omega = omega./sqrt(pi)
xi = sqrt(2)*sigma.*xi + 

% Project 4
% Show monte Carlo and Stochastic collocation together
