function [u, v, p] = LSC_DGS_rotation_2(u, v, p, f1h, f2h, gh, DA1, DA2, B, U, V, width)

%% Initial setup
n = size(u,1); h = width / n; dofu = n*(n+1);


%% Step 1: block Gauss-Seidel relaxation of velocity
[u, v] = symmetricGaussSeidel_Velocity_rotation(u, v, p, f1h, f2h, DA1, DA2, B, U, V);


%% Step 2: Distributive relaxation of velocity and pressue
% rc = gh + div u
vel = [u(:); v(:)];
rc = gh(:) - B*vel;
rc = reshape(rc, n, n);


%% Ap*dq = rc
dq = zeros(n,n); 
dq = symmetricGaussSeidel_Pressure(dq, rc, h);


%% u = u + grad dq
dq = dq(:);
vel = vel + B'*dq;
u = reshape(vel(1:dofu,1), n, n+1);
v = reshape(vel(dofu+1:2*dofu,1), n+1, n);


%%  p = p - inv(Ap) * (-div) * (-mu diff + omega curl) * grad *dq
bigA = [DA1 U; V DA2];
rp = B*bigA*B'*dq;
rp = reshape(rp, n, n);

% Step 4: Ap*(p_{k+1} - p_{k}) = - rp, Ap*(p_{k} - p{k+1}) = rp
ep = zeros(n,n); 
ep = symmetricGaussSeidel_Pressure(ep, rp, h);   
p = p - ep;

end