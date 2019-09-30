function [uh, vh, ph, step, time, residual] = VcycleRotation(uh, vh, ph, f1h, f2h, gh, omega, level, width, nu)

tic

%% Form matrices
A1 = cell(level,1); A2 = cell(level,1); B = cell(level,1);
U = cell(level,1); V = cell(level,1); Lap = cell(level,1);
U2 = cell(level,1); V2 = cell(level,1);
for j = 1:level
    n = 2^(j+1); lef = 0; rig = 1; bot = 0; top = 1; h = (rig - lef) / n; 
    [ux,uy] = meshgrid(lef: h: rig, top-h/2: -h: bot+h/2);
    [vx,vy] = meshgrid(lef+h/2: h: rig-h/2, top: -h: bot);
    omegaU = omega(ux(:),uy(:)); omegaU = reshape(omegaU, n, n+1);
    omegaV = omega(vx(:),vy(:)); omegaV = reshape(omegaV, n+1, n);
    [U{j}, V{j}] = formReactionMatrix(omegaU, omegaV, n);
    
    U2{j} = spdiags(sum(U{j},2), 0, n*(n+1), n*(n+1));
    V2{j} = spdiags(sum(V{j},2), 0, n*(n+1), n*(n+1));
    
    Lap{j} = formDiffusionMatrixPre(width, n);

    A1{j} = formDiffusionMatixU(nu, n, width); 
    A2{j} = formDiffusionMatixV(nu, n, width);
    B{j} = formDivMatrix(n, width);
end
clear ux uy vx vy n lef rig top bot omegaU omegaV


%% Loop
error = 1; step = 0; s1 = 1; s2 = 1; n = size(uh,1); residual = [];
[rh1, rh2, rh3] = getResidualRotation(uh, vh, ph, f1h, f2h, gh, A1{level}, A2{level}, B{level}, U{level}, V{level}); 
while (error > 1e-8)
    euh = zeros(n,n+1); evh = zeros(n+1,n); eph = zeros(n,n);
    [euh, evh, eph] = Vcycle(euh, evh, eph, rh1, rh2, rh3, level);
    uh = uh + euh; vh = vh + evh; ph = ph + eph;
    
    [rh1, rh2, rh3] = getResidualRotation(uh, vh, ph, f1h, f2h, gh, A1{level}, A2{level}, B{level}, U{level}, V{level}); 
    error = sqrt(((norm(rh1(:)))^2 + (norm(rh2(:)))^2 + (norm(rh3(:)))^2)/ ((norm(f1h(:)))^2 + (norm(f2h(:)))^2 + (norm(gh(:)))^2)); 
    step = step + 1; residual(step) = error; 
    fprintf('#Inter step: %8.0u, inter error = %8.4e\n', step, error);
end
ph = ph - mean(ph(:));
time = toc;

    function [u, v, p] = Vcycle(u, v, p, t1, t2, t3, J)
        if J == 1
            [u, v, p] = directSolverOnCoarestMesh(u, v, p, A1{1}, A2{1}, U{1}, V{1}, B{1}, t1, t2, t3);         
            return;
        end
        
        % Pre-smooth
        for i = 1:s1
            [u, v, p] = LSC_DGS_rotation(u, v, p, t1, t2, t3, A1{J}, A2{J}, B{J}, U2{J}, V2{J}, Lap{J}, width);
        end
        
        % Get residual
        [th1, th2, th3] = getResidualRotation(u, v, p, t1, t2, t3, A1{J}, A2{J}, B{J}, U{J}, V{J}); 
        
        % Restrict residual to the coarse grid
        [tH1, tH2, tH3] = Res(th1, th2, th3);
        
        % Cycle
        m = size(tH1,1); duH = zeros(m,m+1); dvH = zeros(m+1,m); dpH = zeros(m,m);
        [duH, dvH, dpH] = Vcycle(duH, dvH, dpH, tH1, tH2, tH3, J-1);
        
        % Prolongate the correction to the fine grid
        [duh, dvh, dph] = Pro(duH, dvH, dpH);
        u = u + duh; v = v + dvh; p = p + dph; 
        
        % Post-smooth
        for i = 1:s2
            [u, v, p] = LSC_DGS_rotation(u, v, p, t1, t2, t3, A1{J}, A2{J}, B{J}, U2{J}, V2{J}, Lap{J}, width);
        end
    end
end