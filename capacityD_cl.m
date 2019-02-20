function[K] = capacityD_cl (x, D, x_cc)
% the function, which calculates maximum possible population of a single
% cluster in the system with a particular phenotype
% x - phenotypic coordinates
% D - dimensionality of the phenotypic space
% x_cc - current location of CCC
S = 0;
for i = 1:D
    S = S + (x(i)-x_cc(i))^4;   
end
K = exp(-S/4);
end