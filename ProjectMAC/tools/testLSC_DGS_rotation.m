function testLSC_DGS_rotation(n, nu)

%% Form matrices
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

clear ux uy vx vy lef rig top bot omegaU omegaV

%% Lopop
error = 1; step = 0; 
[rh1, rh2, rh3] = getResidualRotation(uh, vh, ph, f1h, f2h, gh, A1, A2, B, U, V); 
while (error > 1e-8)
    euh = zeros(n,n+1); evh = zeros(n+1,n); eph = zeros(n,n);
    [euh, evh, eph] = LSC_DGS_rotation(euh, evh, eph, rh1, rh2, rh3, A1, A2, B, U, V, width);
    uh = uh + euh; vh = vh + evh; ph = ph + eph;
    
    [rh1, rh2, rh3] = getResidualRotation(uh, vh, ph, f1h, f2h, gh, A1, A2, B, U, V); 
%     figure(1), surf(rh1); shading interp; colorbar; view(2);
    error = sqrt((norm(rh1(:)))^2 + (norm(rh2(:)))^2 + (norm(rh3(:)))^2)
%     / ((norm(f1h(:)))^2 + (norm(f2h(:)))^2 + (norm(gh(:)))^2))  
    step = step + 1; 
end


end