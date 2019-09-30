function step = drivenCavity(n, J)


%% Setup
h = 2 / n; lef = -1; rig = 1; top = 1; bot = -1;
[ux,uy] = meshgrid(lef: h: rig, top-h/2: -h: bot+h/2);
[vx,vy] = meshgrid(lef+h/2: h: rig-h/2, top: -h: bot);
[px,py] = meshgrid(lef+h/2: h: rig-h/2, top-h/2: -h: bot+h/2);
uh = zeros(n,n+1); vh = zeros(n+1,n); ph = zeros(n,n);

f1h = zeros(n, n+1); f2h = zeros(n+1, n); gh = zeros(n,n);

uTop = ones(1,n+1); uBot = zeros(1,n+1);
vLef = zeros(n+1,1); vRig = zeros(n+1,1);


%% Loop
tol = 1; step = 0; mu = 10;
while (tol > 1e-6)
    [uh, vh, ph] = VCycle(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, J);
    
    [t1, t2, t3] = getResidual(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, h);

    tol = sqrt( (norm(t1,'fro'))^2 + (norm(t2,'fro'))^2 + (norm(t3,'fro'))^2);
    step = step + 1;
end

    function [eu, ev, ep] = VCycle(eu, ev, ep, r1, r2, r3, uT, uB, vL, vR, J)
        if J == 1
            d = size(r1,1);
            eu = zeros(d,d+1); ev = zeros(d+1,d); ep = zeros(d,d); 
            err = 1;
            while (err > 1e-6)
                [eu, ev, ep] = DGS(eu, ev, ep, r1, r2, r3, zeros(1,d+1), zeros(1,d+1), zeros(d+1,1), zeros(d+1,1), 2/d);
                [temp1, temp2, temp3] = getResidual(eu, ev, ep, r1, r2, r3, zeros(1,d+1), zeros(1,d+1), zeros(d+1,1), zeros(d+1,1), 2/d);
                err1 = norm(temp1,'fro'); err4 = norm(r1,'fro');
                err2 = norm(temp2,'fro'); err5 = norm(r2,'fro');
                err3 = norm(temp3,'fro'); err6 = norm(r3,'fro');
                err = sqrt(err1^2 + err2^2 + err3^2) / sqrt(err4^2 + err5^2 + err6^2);
            end
            return;
        end
        
        % Presmooth
        for ii = 1:mu
            [eu, ev, ep] = DGS(eu, ev, ep, r1, r2, r3, uT, uB, vL, vR, 2/size(r1,1));
        end
        
        % Residual
        [R1, R2, R3] = getResidual(eu, ev, ep, r1, r2, r3, uT, uB, vL, vR, 2/size(r1,1));
        
        % Restrict the residual to the coarse grid
        [rH1, rH2, rH3] = Res(R1, R2, R3);
        
        % Vcycle
        d = size(rH1,1);
        duH = zeros(d,d+1); dvH = zeros(d+1,d); dpH = zeros(d,d); 
        [duH, dvH, dpH] = VCycle(duH, dvH, dpH, rH1, rH2, rH3, zeros(1,d+1), zeros(1,d+1), zeros(d+1,1), zeros(d+1,1), J-1);
  
        % Prolongate the correction to the fine grid
        [duh, dvh, dph] = Pro(duH, dvH, dpH);
        eu = eu + duh; ev = ev + dvh; ep = ep + dph;
        
        % Postsmoothing by DGS
        for ii = 1:mu
            [eu, ev, ep] = DGS(eu, ev, ep, r1, r2, r3, uT, uB, vL, vR, 2/size(r1,1));
        end
    end


%% Figure
figure; surf(ux, uy, uh);  colorbar; xlabel('x'); ylabel('y');  shading interp; view(2); title('U'); axis equal; axis off;
figure; surf(vx, vy, vh);  colorbar; xlabel('x'); ylabel('y');  shading interp; view(2); title('V'); axis equal; axis off;
figure; surf(px, py, ph);  colorbar; xlabel('x'); ylabel('y');  shading interp; view(2); title('P'); axis equal; axis off;
% 
figure, contour(flipud(uh),20); xlabel('x'); ylabel('y'); title('U');
figure, contour(flipud(vh),20); xlabel('x'); ylabel('y'); title('V');
figure, contour(flipud(ph),20); xlabel('x'); ylabel('y'); title('P');

end