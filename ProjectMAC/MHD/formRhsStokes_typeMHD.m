function [u1hnew, u2hnew, B1hnew, B2hnew, f1h, f2h, g1h, g2h, dve, dmg] = formRhsStokes_typeMHD(data, ...
                                     lef, rig, bot, top, ux, uy, vx, vy, t, dt, n, h, nu, ...
                                     u1hold, u2hold, B1hold, B2hold)

u1hnew = zeros(n,n+1); u2hnew = zeros(n+1,n); B1hnew = zeros(n,n+1); B2hnew = zeros(n+1,n);
f1h = reshape(data.f1(ux(:),uy(:),t), n, n+1); f2h = reshape(data.f2(vx(:),vy(:),t), n+1, n);
g1h = reshape(data.g1(ux(:),uy(:),t), n, n+1); g2h = reshape(data.g2(vx(:),vy(:),t), n+1, n);

%% Generate boundary data at time = t
xx = lef: h: rig; j = 1:n+1;  yy = top: -h: bot; i = 1:n+1;
u1Lef = data.u1exact(lef,uy(:,1),t); u1Rig = data.u1exact(rig,uy(:,end),t); u1hnew(:,[1 end]) = [u1Lef  u1Rig];
u2Top = data.u2exact(vx(1,:),top,t); u2Bot = data.u2exact(vx(end,:),bot,t); u2hnew([1 end],:) = [u2Top; u2Bot];
u1Top(1,j) = data.u1exact(xx(j),top,t); u1Bot(1,j) = data.u1exact(xx(j),bot,t);
u2Lef(i,1) = data.u2exact(lef,yy(i),t); u2Rig(i,1) = data.u2exact(rig,yy(i),t);

B1Lef = data.B1exact(lef,uy(:,1),t); B1Rig = data.B1exact(rig,uy(:,end),t); B1hnew(:,[1 end]) = [B1Lef  B1Rig];
B2Top = data.B2exact(vx(1,:),top,t); B2Bot = data.B2exact(vx(end,:),bot,t); B2hnew([1 end],:) = [B2Top; B2Bot];
B1Top(1,j) = data.B1exact(xx(j),top,t); B1Bot(1,j) = data.B1exact(xx(j),bot,t);
B2Lef(i,1) = data.B2exact(lef,yy(i),t); B2Rig(i,1) = data.B2exact(rig,yy(i),t);

%% Firstly, add time-dependent term to rhs
% For a boundary node i, mass matrix M(i,i) = 0,
% stiffness matrix A(i,i) = 1, rhs from body force f(i) = gD(x_i,y_i)
% and rhs from nonlinear only take nonzero in interior nodes
f1h = f1h + u1hold/dt; f2h = f2h + u2hold/dt;
g1h = g1h + B1hold/dt; g2h = g2h + B2hold/dt;

%% Secondly, modify rhs for -Lap term
f1h(:,[1 end]) = [u1Lef  u1Rig]; f2h([1 end],:) = [u2Top;  u2Bot]; 
f1h(:,[2 end-1]) = f1h(:,[2 end-1]) + nu/h^2*[u1Lef  u1Rig]; 
f2h([2 end-1],:) = f2h([2 end-1],:) + nu/h^2*[u2Top; u2Bot]; 
j = 2:n; f1h(1,j) = f1h(1,j) + nu/h^2*2*u1Top(1,j); f1h(end,j) = f1h(end,j) + nu/h^2*2*u1Bot(1,j);
i = 2:n; f2h(i,1) = f2h(i,1) + nu/h^2*2*u2Lef(i,1); f2h(i,end) = f2h(i,end) + nu/h^2*2*u2Rig(i,1);

g1h(:,[1 end]) = [B1Lef  B1Rig]; g2h([1 end],:) = [B2Top;  B2Bot]; 
g1h(:,[2 end-1]) = g1h(:,[2 end-1]) + nu/h^2*[B1Lef  B1Rig]; 
g2h([2 end-1],:) = g2h([2 end-1],:) + nu/h^2*[B2Top; B2Bot]; 
j = 2:n; g1h(1,j) = g1h(1,j) + nu/h^2*2*B1Top(1,j); g1h(end,j) = g1h(end,j) + nu/h^2*2*B1Bot(1,j);
i = 2:n; g2h(i,1) = g2h(i,1) + nu/h^2*2*B2Lef(i,1); g2h(i,end) = g2h(i,end) + nu/h^2*2*B2Rig(i,1);

%% Lastly, modify rhs for -div term
% dve is shorthand of divergence of velocity, analogously, dmg
dve = zeros(n,n); dmg = zeros(n,n);
dve(:,[1 end]) = dve(:,[1 end]) + [-u1Lef   u1Rig]/h;  dve([1 end],:) = dve([1 end],:) + [u2Top;  -u2Bot]/h; 
dmg(:,[1 end]) = dmg(:,[1 end]) + [-B1Lef   B1Rig]/h;  dmg([1 end],:) = dmg([1 end],:) + [B2Top;  -B2Bot]/h; 
dve = dve - mean(dve(:));  dmg = dmg - mean(dmg(:));

end