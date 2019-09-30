function [uh, vh, ph] = DefectionCorrectionNS(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, level, au, bu, av, bv, infFlow, nu)

time = 0; t0 = cputime; 

%% Loop
s1 = 1; s2 = 1; alpha_u = 4/3; alpha_p =1; n = size(ph,1); 
[uh, vh, ph] = WCYCLE(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, level);

MaxIt = 6; 
for k = 1:MaxIt
    [r1, r2, r3] = getResidualOseen(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, au{level}, bu{level}, av{level}, bv{level}, 0, nu);
    
    eu = zeros(n,n+1); ev = zeros(n+1,n); ep = zeros(n,n); 
    for j = 1:2
        [eu, ev, ep] = WCYCLE(eu, ev, ep, r1, r2, r3, zeros(1,n+1), zeros(1,n+1), zeros(n+1,1), zeros(n+1,1), level);
    end
  
    uh = uh + eu; vh = vh + ev; ph = ph + ep; 
end
time = time + cputime - t0;


    function [u, v, p] = WCYCLE(u, v, p, f1, f2, f3, uT, uB, vL, vR, J)
        if J == 1
%             err = 1; 
%             while (err > 1e-2)
                [u, v, p] = LSC_DGS(u, v, p, f1, f2, f3, uT, uB, vL, vR, au{1}, bu{1}, av{1}, bv{1}, infFlow, nu);
%                 [temp1, temp2, temp3] = getResidualOseen(u, v, p, f1, f2, f3, uT, uB, vL, vR, au{1}, bu{1}, av{1}, bv{1}, infFlow/(2*size(u,1)));
%                 err1 = (norm(temp1(:)))^2 +  (norm(temp2(:)))^2 + (norm(temp3(:)))^2; err2 = (norm(f1(:)))^2 +  (norm(f2(:)))^2 + (norm(f3(:)))^2;
%                 err = sqrt(err1 / err2);
%             end
            return;
        end
        
        % Relax s1 steps
        for i = 1:s1
            [u, v, p] = LSC_DGS(u, v, p, f1, f2, f3, uT, uB, vL, vR, au{J}, bu{J}, av{J}, bv{J}, infFlow, nu);
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
            [u, v, p] = LSC_DGS(u, v, p, f1, f2, f3, uT, uB, vL, vR, au{J}, bu{J}, av{J}, bv{J}, infFlow, nu);
        end
    end

end