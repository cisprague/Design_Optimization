function [f] = func(x1, x2)
  v1 = x1.^2;
  v2 = x2.*sin(v1);
  f  = v1 + v2;
end

f = func(1, 2)
