function [f1h, f2h, gh] = rhs_remove_Dir_nodes_forStokes(f1h, f2h, gh, uI, vI, pI)
%- To remove the Dirichlet nodes of the rhs.


n = size(f1h, 1);
saveFilename = ['kappa_k_',num2str(n)];
load(saveFilename, 'kappa_k_u', 'kappa_k_v', 'mu', 'gamma', 'uTop', 'uBot', 'vLef', 'vRig')
lef = 0; rig = 1; top = 1; bot = 0; h = (rig - lef) / n;

uI = reshape(uI(:), n, n+1);
vI = reshape(vI(:), n+1, n);
pI = reshape(pI(:), n, n);

uh = zeros(n, n+1);
uh(:,1) = uI(:,1);
uh(:,end) = uI(:,end);

vh = zeros(n+1, n);
vh(1, :) = vI(1, :);
vh(end, :) = vI(end, :);

ph = zeros(n,n);


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
rp(i,j) = - (uh(i,j+1) - uh(i,j))/h - (vh(i,j) - vh(i+1,j))/h;


f1h = f1h - ru;
f2h = f2h - rv;
gh = gh - rp;

end