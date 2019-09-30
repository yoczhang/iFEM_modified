function [velh, magh] = VcycleStokes_typeMHD(level, lef, rig, bot, top, data, n, h, dt, NT, ux, uy, vx, vy)

dofu = n*(n+1); dofp = n^2; width = rig - lef;

%% Form block matrices
A = cell(level,1); B = cell(level,1);  Lap_p = cell(level,1);
for j = 1:level
    m = 2^(j+1); nu = 1;
    [M1, M2] = formMassMatricesUV(m); A1 = formDiffusionMatixU(nu, m, width); A2 = formDiffusionMatixV(nu, m, width);
    A{j} = blkdiag(M1/dt+A1, M2/dt+A2); B{j} = formDivMatrix(m, width); Lap_p{j} = formDiffusionMatrixPre(width, m);
    clear M1 M2 A1 A2 m 
end

%% Form initial u1h, u2h, ph, B1h, B2h, f1h, f2h, g1h, g2h
u1hold = reshape(data.u1exact(ux(:),uy(:),0), n, n+1); u2hold = reshape(data.u2exact(vx(:),vy(:),0), n+1, n);
B1hold = reshape(data.B1exact(ux(:),uy(:),0), n, n+1); B2hold = reshape(data.B2exact(vx(:),vy(:),0), n+1, n);

%% Loop
bigA = [A{level} B{level}'; B{level} sparse(dofp,dofp)]; mu = 1;
fprintf('#Total time step is: %2.0u\n', NT);
for stepTime = 1:NT
    % get iterative initial values by Dirichlet boundary condition at
    % current time and new rhs by 'old' values at last time step
    % attention: nonlinear term is now not included in rhs
    t = stepTime * dt;
    [u1hnew, u2hnew, B1hnew, B2hnew, f1h, f2h, g1h, g2h, dve, dmg] = formRhsStokes_typeMHD(data, ...
                                     lef, rig, bot, top, ux, uy, vx, vy, t, dt, n, h, nu, ...
                                     u1hold, u2hold, B1hold, B2hold);
    phnew = zeros(n,n); rhnew = zeros(n,n);
    i = 1:n; j = 2:n; u1hnew(i,j) = u1hold(i,j); B1hnew(i,j) = B1hold(i,j);
    u2hnew(j,i) = u2hold(j,i); B2hnew(j,i) = B2hold(j,i);
    % get rhs of error equation, so u1,u2,B1,B2 are all homogeneous Dirichlet
    
    
    ERROR = 1; ITE = 0;
    while(ERROR>1e-6)
        [f1, f2, g1, g2] = formNonlinearRhsStokes_typeMHD(data, n, h, u1hnew, u2hnew, B1hnew, B2hnew, t, lef, rig, bot, top);
        bigF = [f1h(:)+f1(:); f2h(:)+f2(:); dve(:)]; bigG = [g1h(:)+g1(:); g2h(:)+g2(:); dmg(:)];
        velh = [u1hnew(:); u2hnew(:); phnew(:)]; magh = [B1hnew(:); B2hnew(:); rhnew(:)];
        
        error = 1; ite = 0;
        while(error>1e-6)
            rhsVelh = bigF - bigA * velh; rhsMagh = bigG - bigA * magh;
            % cycle
            dvelh = zeros(2*dofu+dofp,1); dmagh = zeros(2*dofu+dofp,1);
            [dvelh, dmagh] = Vcycle(dvelh, dmagh, rhsVelh, rhsMagh, level);
            % update
            velh = velh + dvelh; magh = magh + dmagh; ite = ite + 1;
            error = max(norm(bigF - bigA*velh, 'fro') / norm(bigF, 'fro'), norm(bigG - bigA*magh, 'fro') / norm(bigG, 'fro'));
        end
        ITE = ITE + 1;
        u1hnew2 = reshape(velh(1:dofu),n,n+1); u2hnew2 = reshape(velh(dofu+1:2*dofu),n+1,n);
        B1hnew2 = reshape(magh(1:dofu),n,n+1); B2hnew2 = reshape(magh(dofu+1:2*dofu),n+1,n);
        ERROR = max([norm(u1hnew-u1hnew2,'fro')/norm(u1hnew,'fro'), norm(u2hnew-u2hnew2,'fro')^2/norm(u2hnew,'fro'), ...
                     norm(B1hnew-B1hnew2,'fro')/norm(B1hnew,'fro'), norm(B2hnew-B2hnew2,'fro')^2/norm(B2hnew,'fro')]);
        u1hnew = u1hnew2; u2hnew = u2hnew2; B1hnew = B1hnew2; B2hnew = B2hnew2;
    end
    fprintf('#Time step is: %2.0u', stepTime); fprintf(',   Iteration step of cycle is: %2.0u', ite);
    fprintf(',   Iteration step of nonlinear loop is: %2.0u\n', ITE);
    
    u1hold = reshape(velh(1:dofu), n, n+1); u2hold = reshape(velh(dofu+1:2*dofu), n+1, n);
    B1hold = reshape(magh(1:dofu), n, n+1); B2hold = reshape(magh(dofu+1:2*dofu), n+1, n);
end
velh(2*dofu+1:end) = velh(2*dofu+1:end) - mean(velh(2*dofu+1:end)); 
magh(2*dofu+1:end) = magh(2*dofu+1:end) - mean(magh(2*dofu+1:end)); 

    function [du, dB] = Vcycle(du, dB, ru, rB, J)
        if J==1
            leftMatrix = [A{1} B{1}'; B{1} sparse(16,16)]; err = 1;
            while(err>1e-3)
                du = LSC_DGS_Stokes_typeMHD(du, ru, A{1}, B{1}, Lap_p{1});
                dB = LSC_DGS_Stokes_typeMHD(dB, rB, A{1}, B{1}, Lap_p{1});
                err = max(norm(ru-leftMatrix*du,'fro')/norm(ru,'fro'), norm(rB-leftMatrix*dB,'fro')/norm(rB,'fro'));
            end
            return;
        end
        
        % pre-smoothing
        for ii = 1:mu
            du = LSC_DGS_Stokes_typeMHD(du, ru, A{J}, B{J}, Lap_p{J}); dB = LSC_DGS_Stokes_typeMHD(dB, rB, A{J}, B{J}, Lap_p{J});
        end
        
        % restrict residual to coarser mesh
        a = sqrt(size(B{J},1)); dofuh = a*(a+1); dofph = a^2; leftMatrix = [A{J} B{J}'; B{J} sparse(dofph,dofph)];
        ruh = ru - leftMatrix*du; rBh = rB - leftMatrix*dB;
        [ruH1, ruH2, ruH3] = Res(reshape(ruh(1:dofuh),a,a+1), reshape(ruh(dofuh+1:2*dofuh),a+1,a), reshape(ruh(2*dofuh+1:end),a,a));
        [rBH1, rBH2, rBH3] = Res(reshape(rBh(1:dofuh),a,a+1), reshape(rBh(dofuh+1:2*dofuh),a+1,a), reshape(rBh(2*dofuh+1:end),a,a));
        ruH = [ruH1(:);ruH2(:);ruH3(:)]; rBH = [rBH1(:);rBH2(:);rBH3(:)];
        
        % Vcycle
        duH = zeros(a*(a/2+1)+a^2/4,1); dBH = zeros(a*(a/2+1)+a^2/4,1);
        [duH, dBH] = Vcycle(duH, dBH, ruH, rBH, J-1);
        
        % prolongate 
        b = sqrt(size(B{J-1},1)); dof = b*(b+1); 
        [duh1, duh2, duh3] = Pro(reshape(duH(1:dof),b,b+1), reshape(duH(dof+1:2*dof),b+1,b), reshape(duH(2*dof+1:end),b,b));
        [dBh1, dBh2, dBh3] = Pro(reshape(dBH(1:dof),b,b+1), reshape(dBH(dof+1:2*dof),b+1,b), reshape(dBH(2*dof+1:end),b,b));
        duh = [duh1(:);duh2(:);duh3(:)]; dBh = [dBh1(:);dBh2(:);dBh3(:)];
        du = du + duh; dB = dB + dBh;
        
        % post-smoothing
        for ii = 1:mu
            du = LSC_DGS_Stokes_typeMHD(du, ru, A{J}, B{J}, Lap_p{J}); dB = LSC_DGS_Stokes_typeMHD(dB, rB, A{J}, B{J}, Lap_p{J});
        end
    end
end