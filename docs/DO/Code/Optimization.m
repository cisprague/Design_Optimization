clear all; close all;
% Instantiate car with
% nominal values.
TW = Car(              ...
    4.3,               ... r1 m
    0.036,             ... alpha0 rad
    20,                ... Q0
    0.68,              ... omega rad/s
    0.8,               ... r2 m
    0.058,             ... alpha1 rad
    [3, 8].*(2*pi/60), ... omega bounds rad/s
    [0.1, 1.5],        ... r2 bounds m
    [0,0.3],           ... alpha1 bounds rad
    2000,              ... tau
    600,               ... # of samples
    0.05               ... noise
);

[x,fval] = TW.Optimize('particle-swarm')
