function [u, v, p] = LSC_DGS(u, v, p, f1h, f2h, gh, uTop, uBot, vLef, vRig, au, bu, av, bv, Lap, infFlow, nu)
% LSC_DGS(u, v, p, f1, f2, f3, uT, uB, vL, vR, au{J}, bu{J}, av{J}, bv{J}, infFlow, nu);

%% Initial setup
n = size(u,1); h = 1 / n; 


%% Step 1: Gauss-Seidel relaxation of velocity
[u, v] = symmetricGaussSeidel_Velocity(u, v, p, f1h, f2h, uTop, uBot, vLef, vRig, infFlow, nu, au, bu, av, bv);


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


%%  p = p - inv(Ap) * (-div) * (-mu diff + flow div) * grad *dq
% Step 1: \grad dq = [dqu, dqv], completed in last section

% Step 2: (-mu diff + flow div) [dqu, dqv] = [dqu2, dqv2]
nu = (infFlow==0)*nu + (infFlow>0)*h*infFlow/2;

dqu2 = zeros(n,n+1); i = 2:n-1; j = 2:n;
dqu2(i,j) = (nu*(4*dqu(i,j) - dqu(i-1,j) - dqu(i+1,j) - dqu(i,j-1) - dqu(i,j+1)) + 0.5*h*au(i,j).*(dqu(i,j+1) - dqu(i,j-1)) + ...
             0.5*h*bu(i,j).*(dqu(i-1,j) - dqu(i+1,j)) ) / h^2;
         
dqu2(1,j) = (nu*(5*dqu(1,j) - dqu(2,j) - dqu(1,j-1) - dqu(1,j+1)) + 0.5*h*au(1,j).*(dqu(1,j+1) - dqu(1,j-1)) + ...
             0.5*h*bu(1,j).*(-dqu(1,j) - dqu(2,j)) ) / h^2;
         
dqu2(n,j) = (nu*(5*dqu(n,j) - dqu(n-1,j) - dqu(n,j-1) - dqu(n,j+1)) + 0.5*h*au(n,j).*(dqu(n,j+1) - dqu(n,j-1)) + ...
             0.5*h*bu(n,j).*(dqu(n-1,j) + dqu(n,j)) ) / h^2;
 

        
dqv2 = zeros(n+1,n); i = 2:n; j = 2:n-1;
dqv2(i,j) = (nu*(4*dqv(i,j) - dqv(i-1,j) - dqv(i+1,j) - dqv(i,j-1) - dqv(i,j+1)) + 0.5*h*av(i,j).*(dqv(i,j+1) - dqv(i,j-1)) + ...
             0.5*h*bv(i,j).*(dqv(i-1,j) - dqv(i+1,j)) ) / h^2;
         
dqv2(i,1) = (nu*(5*dqv(i,1) - dqv(i-1,1) - dqv(i+1,1) - dqv(i,2)) + 0.5*h*av(i,1).*(dqv(i,1) + dqv(i,2)) + ...
             0.5*h*bv(i,1).*(dqv(i-1,1) - dqv(i+1,1)) ) / h^2;
         
dqv2(i,n) = (nu*(5*dqv(i,n) - dqv(i-1,n) - dqv(i+1,n) - dqv(i,n-1)) + 0.5*h*av(i,n).*(-dqv(i,n) - dqv(i,n-1)) + ...
             0.5*h*bv(i,n).*(dqv(i-1,n) - dqv(i+1,n)) ) / h^2;

% Step 3: (-div) [dqu2, dqv2] = rp
rp = zeros(n,n);  i = 1:n; j = 1:n;
rp(i,j) = - (dqu2(i,j+1) - dqu2(i,j)) / h - (dqv2(i,j) - dqv2(i+1,j)) / h;


% Step 4: Ap*(p_{k+1} - p_{k}) = - rp, Ap*(p_{k} - p{k+1}) = rp
ep = zeros(n,n); 
ep = symmetricGaussSeidel_Pressure(ep, rp, h, Lap);   
p = p - ep;


end