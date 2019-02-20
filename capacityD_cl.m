function[K] = capacityD_cl (x,D, x_cc)
% x_cc - carrying capacity point with D-coordinates
S = 0;
for i = 1:D
    S = S + (x(i)-x_cc(i))^4;   % here is should be c(t)
end
K = exp(-S/4);
end