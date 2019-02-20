# Climate


Functions, required for the simulation:
capacityD_cl - calculates maximum possible population of a single cluster in the system with a particular phenotype.

Parameters:

x - phenotypic coordinates;

D - dimensionality of the phenotypic space;

x_cc - current location of CCC.

competitionD_cl - calculates competition between two clusters x and y.

Parameters:

x,y - phenotypic coordinates of clusters, where x influences on y;

D - dimencionality of the phenotypic space;

b - mutation coefficients;

sigma_comp - width of competition kernel;

x_cc - current location of CCC.

PopulDens_cl - calculates ecological dynamics - changing of clusters density.

Parameters:

t - time period of ecologicaldynamics;

N - population sizes of all clusters;

D - dimensionality of the phenotypic space;

b - mutations coefficients;

m - number of clusters;

X - phenotypic coordinates of clusters, written in a single row (1-D of first cluster, (D+1)-2D of second and so on);

CompetitionMatrix - mxm matrix of competition kernels between all pairs of clusters;

x_cc - current location of CCC.

Phenotypes_cl - calculates evolutionary dynamics, changing of clusters phenotypes.

Parameters:

t - time period of ecologicaldynamics;

X - phenotypic coordinates of clusters, written in a single row (1-D of first cluster, (D+1)-2D of second and so on);

D - dimensionality of the phenotypic space;

b - mutations coefficients;

m - number of clusters;

mD - mxD;

sigma_comp - width of competition kernel;

N - population sizes of all clusters;

x_cc - current location of CCC.

