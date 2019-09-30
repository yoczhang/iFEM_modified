function [u, v, p] = LSC_DGS_rotation(u, v, p, f1h, f2h, gh, A1, A2, B, U, V, Lap, width)

%% Initial setup
n = size(u,1); h = width / n; 


%% Step 1: block Gauss-Seidel relaxation of velocity
[u, v] = symmetricGaussSeidel_Velocity_rotation(u, v, p, f1h, f2h, A1, A2, B, U, V);


%% Step 2: Distributive relaxation of velocity and pressue
% rc = gh + div u
rc = zeros(n, n);  i = 1:n; j = 1:n;
rc(i,j) = gh(i,j) + (u(i,j+1) - u(i,j))/h + (v(i,j) - v(i+1,j))/h;


%% Ap*dq = rc
dq = zeros(n,n); 
dq = symmetricGaussSeidel_Pressure(dq, rc, h, Lap);


%% u = u + grad dq
% \grad dq = [dqu, dqv]
dqu = zeros(n,n+1); i = 1:n; j = 2:n; dqu(i,j) = (dq(i,j) - dq(i,j-1))/h;
dqv = zeros(n+1,n); i = 2:n; j = 1:n; dqv(i,j) = (dq(i-1,j) - dq(i,j))/h;
u = u + dqu; v = v + dqv;


%%  p = p - inv(Ap) * (-div) * (-mu diff + omega curl) * grad *dq
% Step 1: \grad dq = [dqu, dqv], completed in last section

% Step 2: (-mu diff + omega curl) [dqu, dqv] = [dqu2, dqv2]
dquv = [dqu(:); dqv(:)];
matrix = [A1 U; V A2];
dquv2 = matrix*dquv;
dof = n*(n+1);
dqu2 = reshape(dquv2(1:dof,1), n, n+1); dqv2 = reshape(dquv2(dof+1:2*dof,1), n+1, n);

% Step 3: (-div) [dqu2, dqv2] = rp
rp = zeros(n,n);  i = 1:n; j = 1:n;
rp(i,j) = - (dqu2(i,j+1) - dqu2(i,j)) / h - (dqv2(i,j) - dqv2(i+1,j)) / h;


% Step 4: Ap*(p_{k+1} - p_{k}) = - rp, Ap*(p_{k} - p{k+1}) = rp
ep = zeros(n,n); 
ep = symmetricGaussSeidel_Pressure(ep, rp, h, Lap);   
p = p - ep;

end