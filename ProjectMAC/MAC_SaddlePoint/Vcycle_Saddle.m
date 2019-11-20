function [uh, vh, ph, step, time] = Vcycle_Saddle(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, h, width, J, mu, gamma)
 
time = 0; t0 = cputime; n = width / h;

[r1, r2, r3] = getResidual_Saddle(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, h, mu, gamma);

%% Loop
error = 1; step = 0; Nite = 30;
while (error > 1e-6)
    du = zeros(n,n+1); dv = zeros(n+1,n); dp = zeros(n,n);
    [du, dv, dp] = VCYCLE(du, dv, dp, r1, r2, r3, zeros(1,n+1), zeros(1,n+1), zeros(n+1,1), zeros(n+1,1), J);
    
    uh = uh + du; vh = vh + dv; ph = ph + dp;
    [r1, r2, r3] = getResidual_Saddle(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, h, mu, gamma);

    error = sqrt( (norm(r1,'fro'))^2 + (norm(r2,'fro'))^2 + (norm(r3,'fro'))^2);
    step = step + 1;
    
    if mod(step,10) == 0
        disp(['step = ', num2str(step), ' --- ', 'error=', num2str(error)])
    end
end
time = time + cputime - t0;

    function [eu, ev, ep] = VCYCLE(eu, ev, ep, r1, r2, r3, uT, uB, vL, vR, J)
%         disp(['VCYCLE: in the ', num2str(J), '-level'])
        if J == 1
            d = size(r1,1);
            eu = zeros(d,d+1); ev = zeros(d+1,d); ep = zeros(d,d); 
            err = 1;
            flag = 0;
            maxIt = 2000;
            while ((err > 1e-6) && (flag <= maxIt))
                [eu, ev, ep] = DGS_Saddle(eu, ev, ep, r1, r2, r3, zeros(1,d+1), zeros(1,d+1), zeros(d+1,1), zeros(d+1,1), width/d, mu, gamma);
                [temp1, temp2, temp3] = getResidual_Saddle(eu, ev, ep, r1, r2, r3, zeros(1,d+1), zeros(1,d+1), zeros(d+1,1), zeros(d+1,1), width/d, mu, gamma);
                err1 = norm(temp1,'fro') / norm(r1,'fro');
                err2 = norm(temp2,'fro') / norm(r2,'fro');
                err3 = norm(temp3,'fro') / norm(r3,'fro');
                err = sqrt(err1^2 + err2^2 + err3^2);
                flag = flag + 1;
            end
            return;
        end
        
        % Presmooth
        for ii = 1:Nite
            [eu, ev, ep] = DGS_Saddle(eu, ev, ep, r1, r2, r3, uT, uB, vL, vR, width/size(r1,1),mu,gamma);
        end
        
        % Residual
        [R1, R2, R3] = getResidual_Saddle(eu, ev, ep, r1, r2, r3, uT, uB, vL, vR, width/size(r1,1), mu, gamma);
        
        % Restrict the residual to the coarse grid
        [rH1, rH2, rH3] = Res(R1, R2, R3);
        [duH, dvH, dpH] = Res(eu, ev, ep);

        % Vcycle
        d = size(rH1,1);
% %         duH = zeros(d,d+1); dvH = zeros(d+1,d); dpH = zeros(d,d); 
        [duH, dvH, dpH] = VCYCLE(duH, dvH, dpH, rH1, rH2, rH3, zeros(1,d+1), zeros(1,d+1), zeros(d+1,1), zeros(d+1,1), J-1);
        
        % Wcycle
%         [rH12, rH22, rH32] = getResidual_Saddle(duH, dvH, dpH, rH1, rH2, rH3, zeros(1,d+1), zeros(1,d+1), zeros(d+1,1), zeros(d+1,1), width/size(r1,1), mu, gamma);
%         [dduH, ddvH, ddpH] = VCYCLE(zeros(d,d+1), zeros(d+1,d), zeros(d,d), rH12, rH22, rH32, zeros(1,d+1), zeros(1,d+1), zeros(d+1,1), zeros(d+1,1), J-1);
%         duH = duH + dduH; dvH = dvH + ddvH; dpH = dpH + ddpH;
  
        % Prolongate the correction to the fine grid
        [duh, dvh, dph] = Pro(duH, dvH, dpH);
        eu = eu + duh; ev = ev + dvh; ep = ep + dph;
        
        % Postsmoothing by DGS
        for ii = 1:Nite
            [eu, ev, ep] = DGS_Saddle(eu, ev, ep, r1, r2, r3, uT, uB, vL, vR, width/size(r1,1),mu,gamma);
        end
    end

end