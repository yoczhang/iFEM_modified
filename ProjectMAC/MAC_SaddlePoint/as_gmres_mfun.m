function r = as_gmres_mfun(Bh)
%- The input Bh contains f1h, f2h and gh

%- Firstly, we need to get the size of uh, vh and ph
%- f1h: [n x n+1]; f2h: [n+1 x n]; gh: [n x n];
L = length(Bh);
n = (sqrt(3*L+1)-1)/3;

location_f1h = 1 : n*(n+1);
location_f2h = n*(n+1)+1 : 2*n*(n+1);
location_gh = 2*n*(n+1)+1 : L;

f1h = Bh(location_f1h); 
f1h = reshape(f1h, n, n+1);

f2h = Bh(location_f2h);
f2h = reshape(f2h, n+1, n);

gh = Bh(location_gh);
gh = reshape(gh, n, n);

saveFilename = ['kappa_k_',num2str(n)];
load(saveFilename, 'mu', 'gamma')

lef = 0; rig = 1; top = 1; bot = 0; h = (rig - lef) / n; width = rig - lef;
uh = zeros(n,n+1);   vh = zeros(n+1,n); ph = zeros(n,n);
level = log2(n)-1;

[uh, vh, ph, ~, ~] = Vcycle_Saddle(uh, vh, ph, f1h, f2h, gh, zeros(1,n+1), zeros(1,n+1), zeros(n+1,1), zeros(n+1,1), h, width, level, mu, gamma);

r = [uh(:); vh(:); ph(:)];

end 