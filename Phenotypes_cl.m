function dX = Phenotypes_cl (t, X, D, b, m, mD, sigma_comp, N, x_cc) 
% calculates evolutionary dynamics, changing of clusters phenotypes
% D - number of dimensions/phenotype characteristics 
% m - number of clusters, m*D - number of equotions
% x_cc - carrying capacity point with D-coordinates
 
dX = zeros(mD,1);
 
for k = 1:mD                % go through whole vector of equotions
    
    if rem(k,D) == 0;
        r = k/D;
        i = D;
    else
    r = floor (k/D) + 1;    % find index of cluster of this equotion
    i = rem (k,D);          % find index of phenotype coordinate 
    end
    
    %zri = D*(r-1)+i;   % i-th phenotype of r-th cluster - the same with k;
    zr1 = D*(r-1)+1;    % first phenotype of r-th cluster
    rD = r*D;           % end of r-th cluster phenotype variables
    big_thing = 0;
    
     
    for s = 1:m
        
        
        zs1 = D*(s-1)+1;
        zsi = D*(s-1)+i;
        sD = s*D;
        X(zs1:sD);
        X(zr1:rD);
                
        big_thing = big_thing +...
            competitionD_cl(X(zs1:sD), X(zr1:rD), D, b, sigma_comp, x_cc) * N(end, s).*...         % influence on cluster r
            ( sum((b(i,:).').* (X((zs1):sD)-(x_cc.')))  - (X(zsi)-X(k))/sigma_comp(i)^2 - ...
            (X(k)-x_cc(i)).^3  );    
    end
    
    dX(k) = N(end,r) *(1./capacityD_cl(X(zr1:rD), D, x_cc)) * big_thing;
end
    
end
 