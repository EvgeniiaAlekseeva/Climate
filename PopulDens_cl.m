function dN = PopulDens_cl (t, N, D, b, m, X, CompetitionMatrix, x_cc)
% calculates ecological dynamics - changing of clusters density
% X - big vector of phenotype coordinates (1-D of first cluster,
% (D+1)-2D of second and so on);
% m - number of clusters; N - a vector of cluster's density;
% D - number of dimensions; b - matrix of competition coefficients.
% x_cc - carrying capacity point with D-coordinates
 
dN = zeros(m,1);
    for i = 1:m       % i - number of a cluster 
        iD = i*D;
        i1 = D*(i-1) + 1;
        dN(i) = N(i)*(1 - sum(CompetitionMatrix(:,i).*N) / capacityD_cl(X(i1:iD), D, x_cc) ); 
    end
end