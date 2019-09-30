function solu = LSC_DGS_convection(solu, A, B, Lap, t)

%% Initial setup
n = sqrt(size(B,1)); h = 1 / n; dofu = n*(n+1);
p = reshape(solu(2*dofu+1:2*dofu+n^2,1),n,n); 
gh = reshape(t(2*dofu+1:end,1), n, n);

%% Step 1: Gauss-Seidel relaxation of velocity
[u, v] = symmetricGaussSeidel_Velocity(solu, A, B, t, dofu);

%% Step 2: Distributive relaxation of velocity and pressue
% rc = gh + div u
rc = zeros(n, n);  i = 1:n; j = 1:n;
rc(i,j) = gh(i,j) + (u(i,j+1) - u(i,j))/h + (v(i,j) - v(i+1,j))/h;

%% Ap*dq = rc
dq = symmetricGaussSeidel_Pressure(zeros(n,n), rc, h, Lap);

%% u = u + grad dq
% \grad dq = [dqu, dqv]
dqu = zeros(n,n+1); i = 1:n; j = 2:n; dqu(i,j) = (dq(i,j) - dq(i,j-1))/h;
dqv = zeros(n+1,n); i = 2:n; j = 1:n; dqv(i,j) = (dq(i-1,j) - dq(i,j))/h;
u = u + dqu; v = v + dqv;

%%  p = p - inv(Ap) * (-div) * (-mu diff + flow div) * grad *dq
% Step 1: \grad dq = [dqu, dqv], completed in last section

% Step 2: (-mu diff + flow div) [dqu, dqv] = [dqu2, dqv2]
dquv = [dqu(:); dqv(:)]; dquv2 = A*dquv;
dqu2 = reshape(dquv2(1:dofu),n,n+1); dqv2 = reshape(dquv2(dofu+1:2*dofu),n+1,n);

% Step 3: (-div) [dqu2, dqv2] = rp
rp = zeros(n,n);  i = 1:n; j = 1:n;
rp(i,j) = - (dqu2(i,j+1) - dqu2(i,j)) / h - (dqv2(i,j) - dqv2(i+1,j)) / h;

% Step 4: Ap*(p_{k+1} - p_{k}) = - rp, Ap*(p_{k} - p{k+1}) = rp
ep = symmetricGaussSeidel_Pressure(zeros(n,n), rp, h, Lap);   
p = p - ep;

solu = [u(:); v(:); p(:)];
end