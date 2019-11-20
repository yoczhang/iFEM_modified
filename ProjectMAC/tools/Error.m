function [errVelL2, errVelH1, errVelInfi, errPreL2] = Error(uh, vh, ph, uI, vI, pI)

n = size(ph,1);
uh = uh(:); vh = vh(:); ph = ph(:);
uI = uI(:); vI = vI(:); pI = pI(:);
ue = uI - uh; ve = vI - vh; pe = pI - ph;

%% relative L2 error of velocity
% errVelL2 = sqrt((norm(ue)^2 + norm(ve)^2) / (norm(uI)^2 + norm(vI)^2));
errVelL2 = sqrt((1/n*norm(ue))^2 + (1/n*norm(ve))^2);

%% relative H1 error of velocity
nx = n+1;  ny = n;  ex = ones(nx,1);  ey = ones(ny,1);
Tx = spdiags([-ex 2*ex -ex], -1:1, nx, nx);
Ty = spdiags([-ey 2*ey -ey], -1:1, ny, ny);
A = kron(speye(nx),Ty) + kron(Tx,speye(ny));
LEF = 1:ny;   RIG = ny*(nx-1)+1:ny*nx;
TOP = ny+1:ny:ny*(nx-2)+1; BOT = ny*2:ny:ny*(nx-1);
a = false;  a(LEF) = true;  a(RIG) = true; inteNodes = find(~a);
b = zeros(nx*ny,1);  b(TOP) = 1;  b(BOT) = 1;  Ae = spdiags(b, 0, nx*ny, nx*ny);
A = A + Ae;  
temp1 = ue(inteNodes)'*A(inteNodes, inteNodes)*ue(inteNodes);
temp2 = uI(inteNodes)'*A(inteNodes, inteNodes)*uI(inteNodes);
clear nx ny ex ey Tx Ty A a b Ae LEF RIG TOP BOT

nx = n;  ny = n+1;  ex = ones(nx,1);  ey = ones(ny,1);
Tx = spdiags([-ex 2*ex -ex], -1:1, nx, nx);
Ty = spdiags([-ey 2*ey -ey], -1:1, ny, ny);
A = kron(speye(nx),Ty) + kron(Tx,speye(ny));
TOP = 1:ny:ny*(nx-1)+1;  BOT = ny:ny:ny*nx;
LEF = 2:ny-1;  RIG = ny*(nx-1)+2:ny*nx-1;
a = false;  a(TOP) = true;  a(BOT) = true;  inteNodes = find(~a);
b = zeros(nx*ny,1);  b(LEF) = 1;  b(RIG) = 1;  Ae = spdiags(b, 0, nx*ny, nx*ny);
A = A + Ae;  
temp3 = ve(inteNodes)'*A(inteNodes,inteNodes)*ve(inteNodes);
temp4 = vI(inteNodes)'*A(inteNodes,inteNodes)*vI(inteNodes);
clear nx ny ex ey Tx Ty A a b Ae LEF RIG TOP BOT
% errVelH1 = sqrt((temp1 + temp3) / (temp2 + temp4));
errVelH1 = sqrt(temp1 + temp3);

%% Infinity
errVelInfi = max(sqrt(ue.^2 + ve.^2));


%% relative L2 error of pressure
if norm(pI)>0
    errPreL2 = norm(pe) / norm(pI);
else
    errPreL2 = norm(pe);
end

end