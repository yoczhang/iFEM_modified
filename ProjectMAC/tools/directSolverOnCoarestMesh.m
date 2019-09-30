function [u, v, p] = directSolverOnCoarestMesh(u, v, p, A1, A2, U, V, B, t1, t2, t3)
% direct solver on the coarest mesh

%% Inte nodes
n = size(u,1); dofu = n*(n+1); dofp = n^2;

nx = n+1;  ny = n;  LEF = 1:ny;   RIG = ny*(nx-1)+1:ny*nx;
a = false;  a(LEF) = true;  a(RIG) = true; inteNodesU = find(~a);
clear nx ny LEF RIG a

nx = n;  ny = n+1; TOP = 1:ny:ny*(nx-1)+1;  BOT = ny:ny:ny*nx;
a = false;  a(TOP) = true;  a(BOT) = true;  inteNodesV = find(~a);
clear nx ny BOT TOP a

%% 
bigA = [A1  U; V   A2]; 
temp =  sparse(dofp,dofp); 
% temp = spdiags(1e-6*ones(dofp,1), 0, dofp, dofp);
matrix = [bigA B'; B temp]; loadVector = [t1(:);  t2(:);  t3(:)];
solu = [u(:); v(:); p(:)];

loadVector2 = loadVector - matrix*solu;
i = [inteNodesU  inteNodesV+dofu  2*dofu+1:2*dofu+dofp-1];
solu(i) = matrix(i,i) \ loadVector2(i);

u = reshape(solu(1:dofu,1), n, n+1);
v = reshape(solu(dofu+1:2*dofu,1), n+1, n);
p = reshape(solu(2*dofu+1:2*dofu+dofp), n, n);


end