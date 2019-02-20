function dX = Phenotypes_cl (t, X, D, b, m, mD, sigma_comp, N, x_cc) 
% calculates evolutionary dynamics, changing of clusters phenotypes
% t - time period of ecologicaldynamics
% X - phenotypic coordinates of clusters, written in a single row (1-D of first cluster, 
% (D+1)-2D of second and so on);
% D - dimensionality of the phenotypic space
% b - mutations coefficients
% m - number of clusters
% mD - mxD
% sigma_comp - width of competition kernel
% N - population sizes of all clusters
% x_cc - current location of CCC

dX = zeros(mD,1);
 
for k = 1:mD                
    
    if rem(k,D) == 0;
        r = k/D;
        i = D;
    else
    r = floor (k/D) + 1;    
    i = rem (k,D);           
    end
    
    zr1 = D*(r-1)+1;    
    rD = r*D;           
    big_thing = 0;
    
     
    for s = 1:m
        
        
        zs1 = D*(s-1)+1;
        zsi = D*(s-1)+i;
        sD = s*D;
        X(zs1:sD);
        X(zr1:rD);
                
        big_thing = big_thing +...
            competitionD_cl(X(zs1:sD), X(zr1:rD), D, b, sigma_comp, x_cc) * N(end, s).*...         
            ( sum((b(i,:).').* (X((zs1):sD)-(x_cc.')))  - (X(zsi)-X(k))/sigma_comp(i)^2 - ...
            (X(k)-x_cc(i)).^3  );    
    end
    
    dX(k) = N(end,r) *(1./capacityD_cl(X(zr1:rD), D, x_cc)) * big_thing;
end
    
end
 