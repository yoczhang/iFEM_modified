function as_gmres_afunc(U)
%- The input U contains uh, vh and ph

%- Firstly, we need to get the size of uh, vh and ph
%- uh: [n x n+1]; vh: [n+1 x n]; ph: [n x n];
L = length(U);
n = (sqrt(3*L+1)-1)/3;

location_uh = 1 : n*(n+1);
location_vh = n*(n+1)+1 : 2*n*(n+1);
location_ph = 2*n*(n+1)+1 : L;

uh = U(location_uh); 
uh = reshape(uh, n, n+1);

vh = U(location_vh);
vh = reshape(vh, n+1, n);

ph = U(location_ph);
ph = reshape(ph, n, n);




end