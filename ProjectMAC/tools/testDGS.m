function [step, time] = testDGS(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, uI, vI, pI, h)

time = 0; t0 = cputime; 

%% Loop
tol = 1; step = 0;
while (tol > 1e-6)
    uh1 = uh;
    vh1 = vh;
    ph1 = ph;
    [uh, vh, ph] = DGS(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, h);
    [rh1, rh2, rh3] = getResidual(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, h);
    
    tol = sqrt( (norm(rh1,'fro'))^2 + (norm(rh2,'fro'))^2 + (norm(rh3,'fro'))^2);
    step = step + 1;
end
time = time + cputime - t0;

erru = norm(uI - uh(:)) / norm(uI);
errv = norm(vI - vh(:)) / norm(vI);
errp = norm(pI - ph(:)) / norm(pI);


end
