%% Setup
clear all
lef = 0; rig = 1; top = 1; bot = 0; T = 1; level = 4; n = 2^(level+1); 
data = exactSolutionMHD; dofu = n*(n+1); dofp = n^2;
width = rig - lef; h = width/n; NT = n; dt = T/NT;

[ux,uy] = meshgrid(lef: h: rig, top-h/2: -h: bot+h/2); 
[vx,vy] = meshgrid(lef+h/2: h: rig-h/2, top: -h: bot); 
[px,py] = meshgrid(lef+h/2: h: rig-h/2, top-h/2: -h: bot+h/2); 

u1I = data.u1exact(ux(:),uy(:),T); u2I = data.u2exact(vx(:),vy(:),T); 
B1I = data.B1exact(ux(:),uy(:),T); B2I = data.B2exact(vx(:),vy(:),T); 
pI  = data.pexact(px(:),py(:),T);  rI  = zeros(n^2,1);

%% Vcycle
[velh, magh] = VcycleStokes_typeMHD(level, lef, rig, bot, top, ...
                                    data, n, h, dt, NT, ux, uy, vx, vy);

u1h = velh(1:dofu); u2h = velh(dofu+1:2*dofu); ph = velh(2*dofu+1:end);
B1h = magh(1:dofu); B2h = magh(dofu+1:2*dofu); rh = magh(2*dofu+1:end);

velocityL2err = sqrt(h^2*sum((u1I - u1h).^2) + h^2*sum((u2I - u2h).^2));
velocityInf = max(sqrt((u1I - u1h).^2 + (u2I - u2h).^2));
pressureL2err = sqrt(h^2*sum((pI - ph).^2));
magneticL2err = sqrt(h^2*sum((B1I - B1h).^2) + h^2*sum((B2I - B2h).^2));
magneticInf = max(sqrt((B1I - B1h).^2 + (B2I - B2h).^2));
rL2err = sqrt(h^2*sum((rI - rh).^2));

[velocityL2err, velocityInf, pressureL2err;
 magneticL2err, magneticInf, rL2err]

u1h = reshape(u1h, n, n+1); u2h = reshape(u2h, n+1, n);
i = 1:n; j = 1:n; temp1 = zeros(n,n); temp2 = zeros(n,n);
temp1(i,j) = (u1h(i,j+1) - u1h(i,j))/h + (u2h(i,j) - u2h(i+1,j))/h;
B1h = reshape(B1h, n, n+1); B2h = reshape(B2h, n+1, n);
temp2(i,j) = (B1h(i,j+1) - B1h(i,j))/h + (B2h(i,j) - B2h(i+1,j))/h;
[sqrt(h^2*sum((temp1(:)).^2)), sqrt(h^2*sum((temp2(:)).^2))]
