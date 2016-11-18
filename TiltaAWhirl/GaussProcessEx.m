clear all; close all;
% The function
fe = @(x) (x-3).*(x.^3).*((x-6).^4);
% Define sample points. Determine this with LHS
x = [1;2;3.5;5.25;7];
y = fe(x);
% set squared exponential covariance function
covfunc = @covSEiso;
% hyp.cov = [log(1); log(sigma)];
hyp.cov = [0;0];
% set likelihood function to Gaussian
likfunc = @likGauss;
% The noise level. Use sn = 0.05
sn = 1e-16;
hyp.lik = log(sn);
% minimise negative log likelihood function
hyp = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, x, y);
% Use the Guassian process
z = linspace(0,7.5,125)';
[m,s2] = gp(hyp, @infExact, [], covfunc, likfunc, x, y, z);
% plot
f = [m+1.96*sqrt(s2); flipdim(m-1.96*sqrt(s2),1)];
fill([z;flipdim(z,1)], f, [7,7,7]/8);
hold on;
plot(z,m);
plot(x,y,'+');
