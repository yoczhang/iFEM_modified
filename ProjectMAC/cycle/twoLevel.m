function [uh, vh, ph] = twoLevel(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, h)


n = size(uh,1); mu = 10; 

%% Step 1: Presmoothing by DGS
for i = 1:mu
    [uh, vh, ph] = DGS(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, h);
end


%% Step 2: Form residual for momentum and continunity equation

[rh1, rh2, rh3] = getResidual(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, h);


%% Step 3: Restrict the residual to the coarse grid

[rH1, rH2, rH3] = Res(rh1, rh2, rh3);


%% Step 4: Call DGS in the coarse grid till converge
N = n / 2; H = 2*h;
duH = zeros(N,N+1); dvH = zeros(N+1,N); dpH = zeros(N,N);

error = 1;
while (error > 1e-6)
    [duH, dvH, dpH] = DGS(duH, dvH, dpH, rH1, rH2, rH3, zeros(1,N+1), zeros(1,N+1), zeros(N+1,1), zeros(N+1,1), H);
    
    [eH1, eH2, eH3] = getResidual(duH, dvH, dpH, rH1, rH2, rH3, zeros(1,N+1), zeros(1,N+1), zeros(N+1,1), zeros(N+1,1), H);
    
    error = sqrt( (norm(eH1,'fro')^2 + (norm(eH2,'fro'))^2 + (norm(eH3,'fro'))^2)) / ...
            sqrt( (norm(rH1,'fro')^2 + (norm(rH2,'fro'))^2 + (norm(rH3,'fro'))^2));
end


%% Step 5: Prolongate the correction to the fine grid
[duh, dvh, dph] = Pro(duH, dvH, dpH);

uh = uh + duh; vh = vh + dvh; ph = ph + dph;

%% Step 6: Postsmoothing by DGS
for i = 1:mu
    [uh, vh, ph] = DGS(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, h);
end

end