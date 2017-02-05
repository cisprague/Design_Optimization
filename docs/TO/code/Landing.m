clear all; close all;
% Initial state of rocket
S = [100;100;100;5;5;-5;1;0;0;0;0;0;0;10000];
% Target state (soft upright landing, no mass target)
T = [0;0;0;0;0;0;1;0;0;0;0;0;0;NaN];

Falcon = Rocket(S,T,5,1,311,30);

nlp = [];
for i=1:40
    Falcon = Rocket(S,T,5,1,311,i);
    nlp=[nlp;Falcon.NLPDim];
end
plot(nlp)