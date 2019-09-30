function [errVel, errPre] = directSolverForOseen(n, viscosity, infinityNorm, nu)

% if infinityNorm>0 
%     nu = infinityNorm/(2*n); 
% elseif infinityNorm==0 % 
%     nu = viscosity; 
% end
 
[data, uh, vh, ph, f1h, f2h, gh, uLef, uRig, uTop, uBot, vLef, vRig, vTop, vBot, uI, vI, pI] = dataNavierStokes(n, viscosity);
[f1h, f2h] = modifyLoadVectorForLaplaceTerm(f1h, f2h, uLef, uRig, vTop, vBot, uTop, uBot, vLef, vRig, nu);
[au, bu, av, bv] = formConvectionCoefficientMatrices(zeros(n,n+1), zeros(n+1,n), data.flowa, data.flowb, 1);
[f1h, f2h] = modifyLoadVectorForNonlinearTerm(f1h, f2h, uLef, uRig, uTop, uBot, vLef, vRig, vTop, vBot, au, bu, av, bv);


A1 = formDiffusionMatixU(nu, n, 1); A2 = formDiffusionMatixV(nu, n, 1);
B = formDivMatrix(n, 1);
[D1, D2] = formConvectionMatrices(au, bu, av, bv);
A = blkdiag(A1+D1, A2+D2); 

bigA = [A B'; B sparse(n^2,n^2)]; clear A1 A2 D1 D2 B A 
gh = gh - mean(gh(:));
bigF = [f1h(:); f2h(:); gh(:)];
bigu = [uh(:); vh(:); ph(:)];

nx = n+1;  ny = n;  LEF = 1:ny;   RIG = ny*(nx-1)+1:ny*nx;
a = false;  a(LEF) = true;  a(RIG) = true; fixedNodeU = find(a);
clear nx ny LEF RIG a 
nx = n;  ny = n+1; TOP = 1:ny:ny*(nx-1)+1;  BOT = ny:ny:ny*nx;
a = false;  a(TOP) = true;  a(BOT) = true;  fixedNodeV = find(a);
clear nx ny TOP BOT a

fixedNodeP = n^2;
fixedNode = [fixedNodeU'; fixedNodeV'+n*(n+1); fixedNodeP+2*n*(n+1)];
freeNode = setdiff(1:2*n*(n+1)+n^2, fixedNode);
bigF = bigF - bigA*bigu;
bigu(freeNode) = bigA(freeNode,freeNode) \ bigF(freeNode);

uh = reshape(bigu(1:n*(n+1)),n,n+1);
vh = reshape(bigu(n*(n+1)+1:2*n*(n+1)),n+1,n);
ph = reshape(bigu(2*n*(n+1)+1:end),n,n); ph = ph - mean(ph(:));

ue = uh(:) - uI; ve = vh(:) - vI; pe = ph(:) - pI;
uL2 = 1/n*norm(ue);  vL2 = 1/n*norm(ve); 
errVel = sqrt(uL2^2 + vL2^2); errPre = 1/n*norm(pe);

% residual = bigF - bigA*bigu;
% X = ['The viscosity is ', num2str(viscosity)]; disp(X);
% X = ['The errors of velocity and pressure are ', num2str([errVel, errPre])]; disp(X);

uI = reshape(uI, n, n+1); pI = reshape(pI, n, n);
% uResidual = reshape(residual(1:n*(n+1)), n, n+1);
figure, surf(flipud(uh)); view(2); shading interp; colorbar;
figure, surf(flipud(vh)); view(2); shading interp; colorbar;
figure, surf(flipud(ph)); view(2); shading interp; colorbar;
end