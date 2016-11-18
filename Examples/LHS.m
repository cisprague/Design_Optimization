% Example using Latin Hyper Sampling
n = 2; % Design variables
s = 10; % Samples
X = lhsdesign(s,n);

% Bounds
xl = [-1;1];
xu = [2;3.5];

% Scaling
for i = 1:n
    X(:,i) = X(:,i)*(xu(i) - xl(i)) + xl(i);
end
figure;
plot(X(:,1), X(:,2), 'ks')
xlabel('x');
ylabel('y');
title('Latin Hyper Cube Sampling')