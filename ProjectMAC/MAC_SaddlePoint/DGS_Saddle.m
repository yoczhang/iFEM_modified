function [u, v, p] = DGS_Saddle(u, v, p, f1h, f2h, gh, uTop, uBot, vLef, vRig, h, mu, gamma)


%% Initial setup
n = size(u,1);
N = 1;
duTop = zeros(size(uTop));
duBot = zeros(size(uBot));
dvLef = zeros(size(vLef));
dvRig = zeros(size(vRig));

du = zeros(size(u));
dv = zeros(size(v));
%% Step 1: Gauss-Seidel relaxation of velocity
for ite = 1:N
    j = 2:n;
    du(1,j) = ( h^2*f1h(1,j) + mu*(2*uTop(1,j) + u(2,j) + u(1,j-1) + u(1,j+1) - 5*u(1,j)) ...
        - h^2*gamma*u(1,j) - h*(p(1,j)-p(1,j-1)) ...
        + mu*(2*duTop(1,j) + du(2,j) + du(1,j+1) + du(1,j-1) ) ) / (5*mu+h^2*gamma);
    du(n,j) = ( h^2*f1h(n,j) + mu*(2*uBot(1,j) + u(n-1,j) + u(n,j-1) + u(n,j+1) - 5*u(n,j)) ...
        - h^2*gamma*u(n,j) - h*(p(n,j)-p(n,j-1)) ...
        + mu*(2*duBot(1,j) + du(n-1,j) + du(n,j+1) + du(n,j-1) ) ) / (5*mu+h^2*gamma);
    
    i = 2:2:n-1; j = 2:2:n;
    du(i,j) = ( h^2*f1h(i,j) + mu*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - 4*u(i,j)) ...
        - h^2*gamma*u(i,j) - h*(p(i,j)-p(i,j-1)) ...
        + mu*(du(i-1,j) + du(i+1,j) + du(i,j+1) + du(i,j-1) ) ) / (4*mu+h^2*gamma);
    i = 3:2:n-1; j = 3:2:n;
    du(i,j) = ( h^2*f1h(i,j) + mu*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - 4*u(i,j)) ...
        - h^2*gamma*u(i,j) - h*(p(i,j)-p(i,j-1)) ...
        + mu*(du(i-1,j) + du(i+1,j) + du(i,j+1) + du(i,j-1) ) ) / (4*mu+h^2*gamma);
    i = 2:2:n-1; j = 3:2:n;
    du(i,j) = ( h^2*f1h(i,j) + mu*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - 4*u(i,j)) ...
        - h^2*gamma*u(i,j) - h*(p(i,j)-p(i,j-1)) ...
        + mu*(du(i-1,j) + du(i+1,j) + du(i,j+1) + du(i,j-1) ) ) / (4*mu+h^2*gamma);
    i = 3:2:n-1; j = 2:2:n;
    du(i,j) = ( h^2*f1h(i,j) + mu*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - 4*u(i,j)) ...
        - h^2*gamma*u(i,j) - h*(p(i,j)-p(i,j-1)) ...
        + mu*(du(i-1,j) + du(i+1,j) + du(i,j+1) + du(i,j-1) ) ) / (4*mu+h^2*gamma);
    
    i = 2:n;
    dv(i,1) = ( h^2*f2h(i,1) + mu*(2*vLef(i,1) + v(i,2) + v(i-1,1) + v(i+1,1) - 5*v(i,1)) ...
        - h^2*gamma*v(i,1) - h*(p(i-1,1)-p(i,1)) ...
        + mu*(2*dvLef(i,1) + dv(i,2) + dv(i-1,1) + dv(i+1,1)) ) / (5*mu+h^2*gamma);
    dv(i,n) = ( h^2*f2h(i,n) + mu*(2*vRig(i,1) + v(i,n-1) + v(i-1,n) + v(i+1,n) - 5*v(i,n)) ...
        - h^2*gamma*v(i,n) - h*(p(i-1,n)-p(i,n)) ...
        + mu*(2*dvRig(i,1) + dv(i,n-1) + dv(i-1,n) + dv(i+1,n)) ) / (5*mu+h^2*gamma);
    
    i = 2:2:n; j = 2:2:n-1;
    dv(i,j) = ( h^2*f2h(i,j) + mu*(v(i-1,j) + v(i+1,j) + v(i,j-1) + v(i,j+1) - 4*v(i,j)) ...
        - h^2*gamma*v(i,j) - h*(p(i-1,j)-p(i,j)) ...
        + mu*(dv(i-1,j) + dv(i+1,j) + dv(i,j+1) + dv(i,j-1) ) ) / (4*mu+h^2*gamma);
    i = 3:2:n; j = 3:2:n-1;
    dv(i,j) = ( h^2*f2h(i,j) + mu*(v(i-1,j) + v(i+1,j) + v(i,j-1) + v(i,j+1) - 4*v(i,j)) ...
        - h^2*gamma*v(i,j) - h*(p(i-1,j)-p(i,j)) ...
        + mu*(dv(i-1,j) + dv(i+1,j) + dv(i,j+1) + dv(i,j-1) ) ) / (4*mu+h^2*gamma);
    i = 2:2:n; j = 3:2:n-1;
    dv(i,j) = ( h^2*f2h(i,j) + mu*(v(i-1,j) + v(i+1,j) + v(i,j-1) + v(i,j+1) - 4*v(i,j)) ...
        - h^2*gamma*v(i,j) - h*(p(i-1,j)-p(i,j)) ...
        + mu*(dv(i-1,j) + dv(i+1,j) + dv(i,j+1) + dv(i,j-1) ) ) / (4*mu+h^2*gamma);
    i = 3:2:n; j = 2:2:n-1;
    dv(i,j) = ( h^2*f2h(i,j) + mu*(v(i-1,j) + v(i+1,j) + v(i,j-1) + v(i,j+1) - 4*v(i,j)) ...
        - h^2*gamma*v(i,j) - h*(p(i-1,j)-p(i,j)) ...
        + mu*(dv(i-1,j) + dv(i+1,j) + dv(i,j+1) + dv(i,j-1) ) ) / (4*mu+h^2*gamma);
end


%% Step 2: Distributive relaxation of velocity and pressue
% rc = gh + div u
rp = zeros(n, n);
i = 1:n; j = 1:n;
rp(i,j) = gh(i,j) + p(i,j) + (u(i,j+1) - u(i,j))/h + (v(i,j) - v(i+1,j))/h ...
    + (du(i,j+1) - du(i,j))/h + (dv(i,j) - dv(i+1,j))/h;


%% Ap*dp = rp
dp = zeros(n, n); 
for ite = 1:1000
    dp2 = dp;
    
    i = 2:2:n-1; j = 2:2:n-1;
    dp(i,j) = (h^2*rp(i,j) + (1+mu)*(dp(i-1,j) + dp(i+1,j) + dp(i,j-1) + dp(i,j+1)))/(4*(1+mu)+h^2*gamma);
    i = 3:2:n-1; j = 3:2:n-1;
    dp(i,j) = (h^2*rp(i,j) + (1+mu)*(dp(i-1,j) + dp(i+1,j) + dp(i,j-1) + dp(i,j+1)))/(4*(1+mu)+h^2*gamma);
    i = 2:2:n-1; j = 3:2:n-1;
    dp(i,j) = (h^2*rp(i,j) + (1+mu)*(dp(i-1,j) + dp(i+1,j) + dp(i,j-1) + dp(i,j+1)))/(4*(1+mu)+h^2*gamma);
    i = 3:2:n-1; j = 2:2:n-1;
    dp(i,j) = (h^2*rp(i,j) + (1+mu)*(dp(i-1,j) + dp(i+1,j) + dp(i,j-1) + dp(i,j+1)))/(4*(1+mu)+h^2*gamma);
    
    i = 2:n-1;
    dp(i,1) = (h^2*rp(i,1) + (1+mu)*(dp(i-1,1) + dp(i+1,1) + dp(i,2)))/(3*(1+mu)+h^2*gamma); % left 
    dp(i,n) = (h^2*rp(i,n) + (1+mu)*(dp(i-1,n) + dp(i+1,n) + dp(i,n-1)))/(3*(1+mu)+h^2*gamma); % right
    dp(n,i) = (h^2*rp(n,i) + (1+mu)*(dp(n-1,i) + dp(n,i-1) + dp(n,i+1)))/(3*(1+mu)+h^2*gamma); % bottom
    dp(1,i) = (h^2*rp(1,i) + (1+mu)*(dp(2,i) + dp(1,i-1) + dp(1,i+1)))/(3*(1+mu)+h^2*gamma); % top
    
    dp(1,1) = (h^2*rp(1,1) + (1+mu)*(dp(2,1) + dp(1,2)))/(2*(1+mu)+h^2*gamma); % left top
    dp(n,1) = (h^2*rp(n,1) + (1+mu)*(dp(n-1,1) + dp(n,2)))/(2*(1+mu)+h^2*gamma); % left bottom
    dp(1,n) = (h^2*rp(1,n) + (1+mu)*(dp(2,n) + dp(1,n-1)))/(2*(1+mu)+h^2*gamma); % right top
    dp(n,n) = (h^2*rp(n,n) + (1+mu)*(dp(n-1,n) + dp(n,n-1)))/(2*(1+mu)+h^2*gamma); % right bottom
    
    error = max(max(abs(dp - dp2)));
%     error = max(max(abs(dq - dq2)));
%     figure(1); surf(abs(dq - dq2)); colorbar; shading interp;
    if error < 1e-3, break; end
end


%% u = u0 + u + grad dp
i = 1:n; j = 2:n;
u(i,j) = u(i,j) + du(i,j) + (dp(i,j) - dp(i,j-1))/h;

i = 2:n; j = 1:n;
v(i,j) = v(i,j) + dv(i,j) + (dp(i-1,j) - dp(i,j))/h;


%%  p = p - dp - Ap*dp
i = 2:n-1; j = 2:n-1;
p(i,j) = p(i,j) - gamma*dp(i,j) - mu*(4*dp(i,j) - dp(i-1,j) - dp(i+1,j) - dp(i,j-1) - dp(i,j+1)) / h^2;

i = 2:n-1;
p(i,1) = p(i,1) - gamma*dp(i,1) - mu*(3*dp(i,1) - dp(i-1,1) - dp(i+1,1) - dp(i,2)) / h^2; % left
p(i,n) = p(i,n) - gamma*dp(i,n) - mu*(3*dp(i,n) - dp(i-1,n) - dp(i+1,n) - dp(i,n-1)) / h^2; % right
p(n,i) = p(n,i) - gamma*dp(n,i) - mu*(3*dp(n,i) - dp(n-1,i) - dp(n,i-1) - dp(n,i+1)) / h^2; % bottom
p(1,i) = p(1,i) - gamma*dp(1,i) - mu*(3*dp(1,i) - dp(2,i) - dp(1,i-1) - dp(1,i+1)) / h^2; % top

p(1,1) = p(1,1) - gamma*dp(1,1) - mu*(2*dp(1,1) - dp(2,1) - dp(1,2)) / h^2; % left top
p(n,1) = p(n,1) - gamma*dp(n,1) - mu*(2*dp(n,1) - dp(n-1,1) - dp(n,2)) / h^2; % left bottom
p(1,n) = p(1,n) - gamma*dp(1,n) - mu*(2*dp(1,n) - dp(2,n) - dp(1,n-1)) / h^2; % right top
p(n,n) = p(n,n) - gamma*dp(n,n) - mu*(2*dp(n,n) - dp(n-1,n) - dp(n,n-1)) / h^2; % right bottom


end