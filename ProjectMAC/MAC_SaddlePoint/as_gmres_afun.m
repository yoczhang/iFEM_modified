function r = as_gmres_afun(U)
%- The input U contains uh, vh and ph

%- Firstly, we need to get the size of uh, vh and ph
%- uh: [n x n+1]; vh: [n+1 x n]; ph: [n x n];
L = length(U);
n = (sqrt(3*L+1)-1)/3;
lef = 0; rig = 1; top = 1; bot = 0; h = (rig - lef) / n;

saveFilename = ['kappa_k_',num2str(n)];
load(saveFilename, 'kappa_k_u', 'kappa_k_v', 'mu', 'gamma')

location_uh = 1 : n*(n+1);
location_vh = n*(n+1)+1 : 2*n*(n+1);
location_ph = 2*n*(n+1)+1 : L;

uh = U(location_uh); 
uh = reshape(uh, n, n+1);

vh = U(location_vh);
vh = reshape(vh, n+1, n);

ph = U(location_ph);
ph = reshape(ph, n, n);

uTop = uh(1,:);
uBot = uh(end, :);
vLef = vh(:, 1);
vRig = vh(:, end);


%- The following is to compute A*x
ru = zeros(n,n+1); 
rv = zeros(n+1,n); 
rp = zeros(n,n);

i = 2:n-1; 
j = 2:n;
ru(i,j) = (ph(i,j) - ph(i,j-1))/h + gamma*uh(i,j) - mu*(uh(i-1,j) + uh(i+1,j) + uh(i,j-1) + uh(i,j+1) - 4*uh(i,j))/h^2;
ru(1,j) = (ph(1,j) - ph(1,j-1))/h + gamma*uh(1,j) - mu*(2*uTop(1,j) + uh(2,j) + uh(1,j-1) + uh(1,j+1) - 5*uh(1,j))/h^2;
ru(n,j) = (ph(n,j) - ph(n,j-1))/h + gamma*uh(n,j) - mu*(uh(n-1,j) + 2*uBot(1,j) + uh(n,j-1) + uh(n,j+1) - 5*uh(n,j))/h^2;
ru = ru + kappa_k_u.*uh;

i = 2:n; 
j = 2:n-1;
rv(i,j) = (ph(i-1,j) - ph(i,j))/h + gamma*vh(i,j) - mu*(vh(i-1,j) + vh(i+1,j) + vh(i,j-1) + vh(i,j+1) - 4*vh(i,j))/h^2;
rv(i,1) = (ph(i-1,1) - ph(i,1))/h + gamma*vh(i,1) - mu*(vh(i-1,1) + vh(i+1,1) + 2*vLef(i,1) + vh(i,2) - 5*vh(i,1))/h^2;
rv(i,n) = (ph(i-1,n) - ph(i,n))/h + gamma*vh(i,n) - mu*(vh(i-1,n) + vh(i+1,n) + vh(i,n-1) + 2*vRig(i,1) - 5*vh(i,n))/h^2;
rv = rv + kappa_k_v.*vh;
    
i = 1:n; j = 1:n;
rp(i,j) = - ph(i,j) - (uh(i,j+1) - uh(i,j))/h + (vh(i,j) - vh(i+1,j))/h;


%- the final results
r = [ru(:); rv(:); rp(:)];

end