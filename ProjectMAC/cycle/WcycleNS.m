function [uh, vh, ph] = WcycleNS(uh, vh, ph, f1h, f2h, gh, level, au, bu, av, bv, infinityNorm, viscosity, examType)

A1 = cell(level,1); A2  = cell(level,1); A  = cell(level,1); 
OA1 = cell(level,1); OA2  = cell(level,1); OA  = cell(level,1); % for overweighting
D1 = cell(level,1); D2  = cell(level,1); 
B = cell(level,1); Lap = cell(level,1);
for k = 1:level
    n = 2^(k+1); 
%     nu = max(infinityNorm/(2*n), viscosity);
    if infinityNorm>0 
        nu = infinityNorm/(2*n); 
    elseif infinityNorm==0 % 
        nu = viscosity; 
    end
    A1{k} = formDiffusionMatixU(nu, n, 1); A2{k} = formDiffusionMatixV(nu, n, 1);
    OA1{k} = formDiffusionMatixU(nu/2, n, 1); OA2{k} = formDiffusionMatixV(nu/2, n, 1);
    B{k} = formDivMatrix(n, 1);
    [D1{k}, D2{k}] = formConvectionMatrices(au{k}, bu{k}, av{k}, bv{k});
    A{k} = blkdiag(A1{k}+D1{k}, A2{k}+D2{k});
    OA{k} = blkdiag(OA1{k}+D1{k}, OA2{k}+D2{k});
    Lap{k} = formDiffusionMatrixPre(1, n);
end
clear A1 A2 D1 D2 au bu av bv OA1 OA2 

%% Loop
n = size(uh,1); dofu = n*(n+1); dofp = n^2; s1 = 1; s2 = 1; alpha_p =1;
if infinityNorm>0
    alpha_u = 4/3;
elseif infinityNorm==0
    alpha_u = 1;
end
bigA = [A{level} B{level}'; B{level} sparse(dofp,dofp)]; 
bigF = [f1h(:); f2h(:); gh(:)]; solu = [uh(:); vh(:); ph(:)];
t = bigF - bigA*solu;

step = 0; error = 1; 
while (error > 1e-10)
    esolu = zeros(2*dofu+dofp,1); esolu = WCYCLE(esolu, t, level);   
    
    solu = solu + esolu; 
    bigF = [f1h(:); f2h(:); gh(:)]; t = bigF - bigA*solu;
    
    if strcmp(examType,'accurateSolution')
        error = norm(t) / norm(bigF);
    else
        error = norm(t);
    end    
    step = step + 1; 
    fprintf('#Inter step: %8.0u, inter error = %8.4e\n', step, error);
end
uh = reshape(solu(1:dofu), n, n+1);
vh = reshape(solu(dofu+1:2*dofu), n+1, n);
ph = reshape(solu(2*dofu+1:end), n, n); ph = ph - mean(ph(:));

    function u = WCYCLE(u, r, J)
        if J == 1
                u = LSC_DGS_convection(u, A{1}, B{1}, Lap{1}, r);
            return;
        end
        
        % Relax s1 steps
        for i = 1:s1
            u = LSC_DGS_convection(u, A{J}, B{J}, Lap{J}, r);
        end
        
        % Get residual
        m = sqrt(size(B{J},1)); bigA1 = [A{J} B{J}'; B{J} sparse(m^2,m^2)];
        rh = r - bigA1*u; clear bigA1
        
        % Restrict the residual to the coarse grid
        rh1 = reshape(rh(1:m*(m+1)),m,m+1);
        rh2 = reshape(rh(m*(m+1)+1:2*m*(m+1)),m+1,m);
        rh3 = reshape(rh(2*m*(m+1)+1:end),m,m);
        [rH1, rH2, rH3] = Res(rh1, rh2, rh3);
        rH = [rH1(:); rH2(:); rH3(:)];
        clear rh1 rh2 rh3 rH1 rH2 rH3
        
        % Wcycle
        z = m/2; du = zeros(2*z*(z+1)+z^2,1);
        du = WCYCLE(du, rH, J-1);
        
        % Modify the right hand side of the coarse grid problem
        if infinityNorm>0 % do overweighting
            bigOA = [A{J-1}-alpha_u*OA{J-1} B{J-1}'-alpha_u*B{J-1}'; 
                     B{J-1}-alpha_p*B{J-1}  sparse(z^2,z^2)];
            rH = rH + bigOA*du;
            clear bigOA
        end
        
        % Wcycle
        bigA2 = [A{J-1} B{J-1}'; B{J-1} sparse(z^2,z^2)];
        rH2 = rH - bigA2*du; clear bigA2
        du2 = zeros(2*z*(z+1)+z^2,1); 
        du2 = WCYCLE(du2, rH2, J-1);
        
        % Prolongate the correction to the fine grid
        du2 = du2 + du; 
        duH2 = reshape(du2(1:z*(z+1)),z,z+1);
        dvH2 = reshape(du2(z*(z+1)+1:2*z*(z+1)),z+1,z);
        dpH2 = reshape(du2(2*z*(z+1)+1:end),z,z);
        [duh2, dvh2, dph2] = Pro(duH2, dvH2, dpH2);
        duh = [duh2(:); dvh2(:); dph2(:)];
        temp = spdiags([alpha_u*ones(2*m*(m+1),1); alpha_p*ones(m^2,1)], 0, 2*m*(m+1)+m^2, 2*m*(m+1)+m^2);
        u = u + temp*duh;  
        clear du2 duH2 dvH2 dpH2 duh2 dvh2 dph2 duh temp

        % Post-smoothing
        for i = 1:s2
            u = LSC_DGS_convection(u, A{J}, B{J}, Lap{J}, r);
        end
    end

end