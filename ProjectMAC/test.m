function test(n, nu)
% [lamdaLarge, lamdaSmall] = 
% some tests for new ideas and basic iteration

%% Using new index system for variable v, newIndexV
% order: 12 11 10; 9 8 7; 6 5 4; 3 2 1

% old = 1:n*(n+1); 
% temp = false(1,n*(n+1)); temp(mod(old,n)==0) = true; temp2 = find(~temp);
% colume = n + 1 - mod(temp2,n); row = n + 1 - floor(temp2/n);
% newIndexV = zeros(1,n*(n+1)); temp3 = row + (colume - 1)*(n + 1);
% newIndexV(temp3) = temp2; newIndexV(1,1:n+1) = n*(n+1): -n: n;
% clear old temp temp2 colume row temp3

% i = 1:n*(n+1); U2(i, newIndexV) = U; V2(newIndexV, i) = V;

%% Form matrices, A1, U, V, A2, B, A, Lap
[data, uh, vh, ph, f1h, f2h, gh, uI, vI, pI, width, nu] = dataOseenRotation(n, nu);

lef = 0; rig = 1; bot = 0; top = 1; h = (rig - lef) / n; 
[ux,uy] = meshgrid(lef: h: rig, top-h/2: -h: bot+h/2);
[vx,vy] = meshgrid(lef+h/2: h: rig-h/2, top: -h: bot);
omegaU = data.omega(ux(:),uy(:)); omegaU = reshape(omegaU, n, n+1);
omegaV = data.omega(vx(:),vy(:)); omegaV = reshape(omegaV, n+1, n);
[U, V] = formReactionMatrix(omegaU, omegaV, n);

A1 = formDiffusionMatixU(nu, n, width); 
A2 = formDiffusionMatixV(nu, n, width);
B = formDivMatrix(n, width);
Lap = formDiffusionMatrixPre(width, n);
A = [A1 U; V A2];
clear ux uy vx vy omegaU omegaV 
    

%% Interior nodes for u and v, inteNodesU, inteNodesV
dofu = n*(n+1); dofp = n^2;
nx = n+1;  ny = n;  LEF = 1:ny;   RIG = ny*(nx-1)+1:ny*nx;
a = false;  a(LEF) = true;  a(RIG) = true; inteNodesU = find(~a);
clear nx ny LEF RIG a

nx = n;  ny = n+1; TOP = 1:ny:ny*(nx-1)+1;  BOT = ny:ny:ny*nx;
a = false;  a(TOP) = true;  a(BOT) = true;  inteNodesV = find(~a);
clear nx ny BOT TOP a


%% Eigenvalue of A, lamdaLarge, lamdaSmall
% freeDof = [inteNodesU inteNodesV+dofu];
% A2 = A(freeDof, freeDof); k = 4;
% lamdaLarge = eigs(A2, k, 'lm');
% lamdaSmall = eigs(A2, k, 'sm');
% % figure, plot(real(lamda), imag(lamda), 'ro', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');


%% LSC-DGS smoother
error = 1; ite = 0;
while(error > 1e-6)
    [uh, vh, ph] = LSC_DGS_rotation(uh, vh, ph, f1h, f2h, gh, A1, A2, B, U, V, Lap, width);

    [rh1, rh2, rh3] = getResidualRotation(uh, vh, ph, f1h, f2h, gh, A1, A2, B, U, V); 
    error = sqrt(((norm(rh1(:)))^2 + (norm(rh2(:)))^2 + (norm(rh3(:)))^2)/ ((norm(f1h(:)))^2 + (norm(f2h(:)))^2 + (norm(gh(:)))^2));  
    ite = ite + 1;
    X2 = ['error = ', num2str(sqrt(((norm(rh1(:)))^2 + (norm(rh2(:)))^2)/ ((norm(f1h(:)))^2 + (norm(f2h(:)))^2)))]; disp(X2);
end
X = ['The viscosity is ', num2str(nu)]; disp(X);
X = ['The iteration step is ', num2str(ite)]; disp(X);


%% Direct solver
% bigA = [A B'; B sparse(dofp,dofp)];
% bigF = [f1h(:); f2h(:); gh(:)];
% solu = [uh(:); vh(:); ph(:)];
% 
% i = [inteNodesU  inteNodesV+dofu  2*dofu+1:2*dofu+dofp-1];
% solu(i) = bigA(i,i) \ bigF(i);
% uh = reshape(solu(1:dofu), n, n+1);
% vh = reshape(solu(dofu+1:2*dofu), n+1, n);
% ph = reshape(solu(2*dofu+1:2*dofu+dofp), n, n);


%% Uzawa iteration
% error = 1; ite = 0; tau = 0.8;
% f = [f1h(:); f2h(:)];
% while(error > 1e-6)
%     ph = ph(:); 
%     bigu = A \ (f - B'*ph);
%     ph = ph - tau * (gh(:) - B*bigu);
%     
%     uh = reshape(bigu(1:dofu,1), n, n+1);
%     vh = reshape(bigu(dofu+1:2*dofu,1), n+1, n); 
%     ph = reshape(ph, n, n);
%     [rh1, rh2, rh3] = getResidualRotation(uh, vh, ph, f1h, f2h, gh, A1, A2, B, U, V); 
%     error = sqrt(((norm(rh1(:)))^2 + (norm(rh2(:)))^2 + (norm(rh3(:)))^2)/ ((norm(f1h(:)))^2 + (norm(f2h(:)))^2 + (norm(gh(:)))^2));  
%     ite = ite + 1;
%     X2 = ['error = ', num2str(sqrt(((norm(rh1(:)))^2 + (norm(rh2(:)))^2)/ ((norm(f1h(:)))^2 + (norm(f2h(:)))^2)))]; disp(X2);
% end
% X = ['The viscosity is ', num2str(nu)]; disp(X);
% X = ['The iteration step is ', num2str(ite)]; disp(X);

%% Figure
% ph = ph - mean(ph(:));
% 
% uL2 = norm(uh(:) - uI) / norm(uI);
% vL2 = norm(vh(:) - vI) / norm(vI);
% pL2 = norm(ph(:) - pI) / norm(pI);
% [uL2, vL2, pL2]

figure, surf(flipud(uh)); shading interp; colorbar; view(2);
figure, surf(flipud(vh)); shading interp; colorbar; view(2);
figure, surf(flipud(ph)); shading interp; colorbar; view(2);
end