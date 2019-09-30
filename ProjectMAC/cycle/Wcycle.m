function [uh, vh, ph, step, time, residual] = Wcycle(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, level, flowa, flowb, infFlow, nu, examType)

time = 0; t0 = cputime; 

au = cell(level,1); av = cell(level,1); bu = cell(level,1); bv = cell(level,1); Lap = cell(level,1);
for j = 1:level
    n = 2^(j+1); h = 1 / n; 
    [ux,uy] = meshgrid(0: h: 1, 1-h/2: -h: h/2);
    [vx,vy] = meshgrid(h/2: h: 1-h/2, 1: -h: 0);
    au{j} = flowa(ux(:), uy(:)); bu{j} = flowb(ux(:), uy(:));
    av{j} = flowa(vx(:), vy(:)); bv{j} = flowb(vx(:), vy(:));
    au{j} = reshape(au{j},n,n+1); bu{j} = reshape(bu{j},n,n+1);
    av{j} = reshape(av{j},n+1,n); bv{j} = reshape(bv{j},n+1,n);
    Lap{j} = formDiffusionMatrixPre(1, n);
end
clear ux uy vx vy

%% Loop
n = size(uh,1); s1 = 1; s2 = 1; alpha_u = 4/3; alpha_p =1; 
[t1, t2, t3] = getResidualOseen(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, au{level}, bu{level}, av{level}, bv{level}, infFlow, nu);

step = 0; error = 1; residual = [];
while (error > 1e-6)
    eu = zeros(n,n+1); ev = zeros(n+1,n); ep = zeros(n,n);
    [eu, ev, ep] = WCYCLE(eu, ev, ep, t1, t2, t3, zeros(1,n+1), zeros(1,n+1), zeros(n+1,1), zeros(n+1,1), level);   
    
    uh = uh + eu; vh = vh + ev; ph = ph + ep;
     
    [t1, t2, t3] = getResidualOseen(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, au{level}, bu{level}, av{level}, bv{level}, infFlow, nu);
    
    if strcmp(examType,'accurateSolution')
        error = sqrt(((norm(t1(:)))^2 + (norm(t2(:)))^2 + (norm(t3(:)))^2) / ((norm(f1h(:)))^2 + (norm(f2h(:)))^2 + (norm(gh(:)))^2));
    else
        error = sqrt((norm(t1(:)))^2 + (norm(t2(:)))^2 + (norm(t3(:)))^2);
    end
    error
    step = step + 1; residual(step) = error;
end
time = time + cputime - t0;


    function [u, v, p] = WCYCLE(u, v, p, f1, f2, f3, uT, uB, vL, vR, J)
        if J == 1
%             err = 1;
%             while (err > 1e-2)
                [u, v, p] = LSC_DGS(u, v, p, f1, f2, f3, uT, uB, vL, vR, au{1}, bu{1}, av{1}, bv{1}, Lap{J}, infFlow, nu);
%                 [temp1, temp2, temp3] = getResidualOseen(u, v, p, f1, f2, f3, uT, uB, vL, vR, au{1}, bu{1}, av{1}, bv{1}, infFlow/(2*size(u,1)));
%                 err1 = (norm(temp1(:)))^2 +  (norm(temp2(:)))^2 + (norm(temp3(:)))^2; err2 = (norm(f1(:)))^2 +  (norm(f2(:)))^2 + (norm(f3(:)))^2;
%                 err = sqrt(err1 / err2);
%             end
            return;
        end
        
        % Relax s1 steps
        for i = 1:s1
            [u, v, p] = LSC_DGS(u, v, p, f1, f2, f3, uT, uB, vL, vR, au{J}, bu{J}, av{J}, bv{J}, Lap{J}, infFlow, nu);
        end
        
        % Get residual
        [rh1, rh2, rh3] = getResidualOseen(u, v, p, f1, f2, f3, uT, uB, vL, vR, au{J}, bu{J}, av{J}, bv{J}, infFlow, nu);
        
        % Restrict the residual to the coarse grid
        [rH1, rH2, rH3] = Res(rh1, rh2, rh3);
        
        % Wcycle
        m = size(rH1,1); du = zeros(m,m+1); dv = zeros(m+1,m); dp = zeros(m,m);
        [du, dv, dp] = WCYCLE(du, dv, dp, rH1, rH2, rH3, zeros(1,m+1), zeros(1,m+1), zeros(m+1,1), zeros(m+1,1), J-1);
        
        % Modify the right hand side of the coarse grid problem
        [rH4, rH5, rH6] = getNewRHS(du, dv, dp, alpha_u, alpha_p, au{J-1}, bu{J-1}, av{J-1}, bv{J-1}, infFlow);
        rH7 = rH4 + rH1; rH8 = rH5 + rH2; rH9 = rH6 + rH3;
        
        % Wcycle
        mm = size(rH4,1);
        [rHu, rHv, rHp] = getResidualOseen(du, dv, dp, rH7, rH8, rH9, zeros(1,mm+1), zeros(1,mm+1), zeros(mm+1,1), zeros(mm+1,1), ...
                                           au{J-1}, bu{J-1}, av{J-1}, bv{J-1}, infFlow, nu);
        duH = zeros(mm,mm+1); dvH = zeros(mm+1,mm); dpH = zeros(mm,mm);
        [duH, dvH, dpH] = WCYCLE(duH, dvH, dpH, rHu, rHv, rHp, zeros(1,mm+1), zeros(1,mm+1), zeros(mm+1,1), zeros(mm+1,1), J-1);
        
        % Prolongate the correction to the fine grid
        duH = duH + du; dvH = dvH + dv; dpH = dpH + dp;
        [duh, dvh, dph] = Pro(duH, dvH, dpH);
        u = u + alpha_u*duh; v = v + alpha_u*dvh; p = p + alpha_p*dph;     

        % Post-smoothing
        for i = 1:s2
            [u, v, p] = LSC_DGS(u, v, p, f1, f2, f3, uT, uB, vL, vR, au{J}, bu{J}, av{J}, bv{J}, Lap{J}, infFlow, nu);
        end
    end

end