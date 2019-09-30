function UzawaRotationForm(n, tau, nu)

%% Form matrices
[data, uh, vh, ph, f1h, f2h, gh, uI, vI, pI, width, nu] = dataOseenRotation(n, nu);

lef = 0; rig = 1; bot = 0; top = 1; h = (rig - lef) / n; dofu = n*(n+1); dofp = n^2;
[ux,uy] = meshgrid(lef: h: rig, top-h/2: -h: bot+h/2);
[vx,vy] = meshgrid(lef+h/2: h: rig-h/2, top: -h: bot);
omegaU = data.omega(ux(:),uy(:)); omegaU = reshape(omegaU, n, n+1);
omegaV = data.omega(vx(:),vy(:)); omegaV = reshape(omegaV, n+1, n);
[U, V] = formReactionMatrix(omegaU, omegaV, n);

A1 = formDiffusionMatixU(nu, n, width); 
A2 = formDiffusionMatixV(nu, n, width);
B = formDivMatrix(n, width);

bigA = [A1 U; V A2]; bigF = [f1h(:); f2h(:)]; gh = gh(:);


%% Uzawa
error = 1; 
while(error > 1e-6)
    ph = ph(:); 
    bigu = bigA \ (bigF - B'*ph);
    ph = ph - tau * (gh - B*bigu);
    
    uh = reshape(bigu(1:dofu,1), n, n+1);
    vh = reshape(bigu(dofu+1:2*dofu,1), n+1, n); 
    ph = reshape(ph, n, n);
    [rh1, rh2, rh3] = getResidualRotation(uh, vh, ph, f1h, f2h, gh, A1, A2, B, U, V); 
    error = sqrt(((norm(rh1(:)))^2 + (norm(rh2(:)))^2 + (norm(rh3(:)))^2)/ ((norm(f1h(:)))^2 + (norm(f2h(:)))^2 + (norm(gh(:)))^2));  
    
    [sqrt(((norm(rh1(:)))^2 + (norm(rh2(:)))^2)/ ((norm(f1h(:)))^2 + (norm(f2h(:)))^2)), norm(rh3(:))]
end


%% Direct solver
% leftMatrix = [bigA B'; B sparse(dofp,dofp)];
% loadVector = [bigf; gh];
% 
% r = leftMatrix \ loadVector;
% 
% uh = reshape(r(1:dofu,1), n, n+1);
% vh = reshape(r(dofu+1:2*dofu,1), n+1, n);
% ph = reshape(r(2*dofu+1:2*dofu+dofp,1), n, n);


%% Figure
% figure, surf(flipud(uh)); shading interp; colorbar; view(2);
% figure, surf(flipud(vh)); shading interp; colorbar; view(2);
% figure, surf(flipud(ph)); shading interp; colorbar; view(2);
end