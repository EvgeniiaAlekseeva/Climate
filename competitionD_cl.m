function[A] = competitionD_cl(x, y, D, b, sigma_comp, x_cc)   
% the function calculates competition between two clusters
% x,y - phenotypic coordinates of clusters, where x influences on y
% D - dimencionality of the phenotypic space
% b - mutation coefficients
% sigma_comp - width of competition kernel
% x_cc - current location of CCC
 
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
