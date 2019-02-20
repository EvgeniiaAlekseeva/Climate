function[A] = competitionD_cl(x,y,D,b,sigma_comp, x_cc)   
% the function calculates competition between two D-dimension phenotype clusters
% x,y - clusters, where x influences on y
% D - number of dimensions
% b - mutation coefficients
% x_cc - carrying capacity point with D-coordinates
 
B = 0;
S = 0;
for i = 1:D
    for j = 1:D
        B = B + b(i,j)*(x(i)-y(i))*(x(j)-x_cc(j));
    end
    S = S + (x(i)-y(i))^2/(2*sigma_comp(i)^2);
end
A = exp(B - S);
end
