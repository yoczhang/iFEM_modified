function [u, v, p] = DGS(u, v, p, f1h, f2h, gh, uTop, uBot, vLef, vRig, h)


%% Initial setup
n = size(u,1); mu = 1; 
  
%% Step 1: Gauss-Seidel relaxation of velocity
for ite = 1:mu
    j = 2:n;
    u(1,j) = (2*uTop(1,j) + u(2,j) + u(1,j-1) + u(1,j+1) - h*(p(1,j) - p(1,j-1)) + h^2*f1h(1,j))/5;
    u(n,j) = (u(n-1,j) + 2*uBot(1,j) + u(n,j-1) + u(n,j+1) - h*(p(n,j) - p(n,j-1)) + h^2*f1h(n,j))/5;
    
    i = 2:2:n-1; j = 2:2:n;
    u(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - h*(p(i,j) - p(i,j-1)) + h^2*f1h(i,j))/4;
    i = 3:2:n-1; j = 3:2:n;
    u(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - h*(p(i,j) - p(i,j-1)) + h^2*f1h(i,j))/4;
    i = 2:2:n-1; j = 3:2:n;
    u(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - h*(p(i,j) - p(i,j-1)) + h^2*f1h(i,j))/4;
    i = 3:2:n-1; j = 2:2:n;
    u(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - h*(p(i,j) - p(i,j-1)) + h^2*f1h(i,j))/4;
    
    i = 2:n;
    v(i,1) = (v(i-1,1) + v(i+1,1) + 2*vLef(i,1) + v(i,2) - h*(p(i-1,1) - p(i,1)) + h^2*f2h(i,1))/5;
    v(i,n) = (v(i-1,n) + v(i+1,n) + v(i,n-1) + 2*vRig(i,1) - h*(p(i-1,n) - p(i,n)) + h^2*f2h(i,n))/5;
    
    i = 2:2:n; j = 2:2:n-1;
    v(i,j) = (v(i-1,j) + v(i+1,j) + v(i,j-1) + v(i,j+1) - h*(p(i-1,j) - p(i,j)) + h^2*f2h(i,j))/4;
    i = 3:2:n; j = 3:2:n-1;
    v(i,j) = (v(i-1,j) + v(i+1,j) + v(i,j-1) + v(i,j+1) - h*(p(i-1,j) - p(i,j)) + h^2*f2h(i,j))/4;
    i = 2:2:n; j = 3:2:n-1;
    v(i,j) = (v(i-1,j) + v(i+1,j) + v(i,j-1) + v(i,j+1) - h*(p(i-1,j) - p(i,j)) + h^2*f2h(i,j))/4;
    i = 3:2:n; j = 2:2:n-1;
    v(i,j) = (v(i-1,j) + v(i+1,j) + v(i,j-1) + v(i,j+1) - h*(p(i-1,j) - p(i,j)) + h^2*f2h(i,j))/4;
end


%% Step 2: Distributive relaxation of velocity and pressue
% rc = gh + div u
rc = zeros(n, n);
i = 1:n; j = 1:n;
rc(i,j) = gh(i,j) + (u(i,j+1) - u(i,j))/h + (v(i,j) - v(i+1,j))/h;


%% Ap*dq = rc
dq = zeros(n, n); 
for ite = 1:1000
    dq2 = dq;
    
    i = 2:2:n-1; j = 2:2:n-1;
    dq(i,j) = (h^2*rc(i,j) + dq(i-1,j) + dq(i+1,j) + dq(i,j-1) + dq(i,j+1))/4;
    i = 3:2:n-1; j = 3:2:n-1;
    dq(i,j) = (h^2*rc(i,j) + dq(i-1,j) + dq(i+1,j) + dq(i,j-1) + dq(i,j+1))/4;
    i = 2:2:n-1; j = 3:2:n-1;
    dq(i,j) = (h^2*rc(i,j) + dq(i-1,j) + dq(i+1,j) + dq(i,j-1) + dq(i,j+1))/4;
    i = 3:2:n-1; j = 2:2:n-1;
    dq(i,j) = (h^2*rc(i,j) + dq(i-1,j) + dq(i+1,j) + dq(i,j-1) + dq(i,j+1))/4;
    
    i = 2:n-1;
    dq(i,1) = (dq(i-1,1) + dq(i+1,1) + dq(i,2) + h^2*rc(i,1))/3; % left 
    dq(i,n) = (dq(i-1,n) + dq(i+1,n) + dq(i,n-1) + h^2*rc(i,n))/3; % right
    dq(n,i) = (dq(n-1,i) + dq(n,i-1) + dq(n,i+1) + h^2*rc(n,i))/3; % bottom
    dq(1,i) = (dq(2,i) + dq(1,i-1) + dq(1,i+1) + h^2*rc(1,i))/3; % top
    
    dq(1,1) = (dq(2,1) + dq(1,2) + h^2*rc(1,1))/2; % left top
    dq(n,1) = (dq(n-1,1) + dq(n,2) + h^2*rc(n,1))/2; % left bottom
    dq(1,n) = (dq(2,n) + dq(1,n-1) + h^2*rc(1,n))/2; % right top
    dq(n,n) = (dq(n-1,n) + dq(n,n-1) + h^2*rc(n,n))/2; % right bottom
    
    error = max(max(abs(dq - dq2)));
%     error = max(max(abs(dq - dq2)));
%     figure(1); surf(abs(dq - dq2)); colorbar; shading interp;
    if error < 1e-3, break; end
end


%% u = u + grad dq
i = 1:n; j = 2:n;
u(i,j) = u(i,j) + (dq(i,j) - dq(i,j-1))/h;

i = 2:n; j = 1:n;
v(i,j) = v(i,j) + (dq(i-1,j) - dq(i,j))/h;


%%  p = p - Ap*dq
i = 2:n-1; j = 2:n-1;
p(i,j) = p(i,j) - (4*dq(i,j) - dq(i-1,j) - dq(i+1,j) - dq(i,j-1) - dq(i,j+1)) / h^2;

i = 2:n-1;
p(i,1) = p(i,1) - (3*dq(i,1) - dq(i-1,1) - dq(i+1,1) - dq(i,2)) / h^2; % left
p(i,n) = p(i,n) - (3*dq(i,n) - dq(i-1,n) - dq(i+1,n) - dq(i,n-1)) / h^2; % right
p(n,i) = p(n,i) - (3*dq(n,i) - dq(n-1,i) - dq(n,i-1) - dq(n,i+1)) / h^2; % bottom
p(1,i) = p(1,i) - (3*dq(1,i) - dq(2,i) - dq(1,i-1) - dq(1,i+1)) / h^2; % top

p(1,1) = p(1,1) - (2*dq(1,1) - dq(2,1) - dq(1,2)) / h^2; % left top
p(n,1) = p(n,1) - (2*dq(n,1) - dq(n-1,1) - dq(n,2)) / h^2; % left bottom
p(1,n) = p(1,n) - (2*dq(1,n) - dq(2,n) - dq(1,n-1)) / h^2; % right top
p(n,n) = p(n,n) - (2*dq(n,n) - dq(n-1,n) - dq(n,n-1)) / h^2; % right bottom


end