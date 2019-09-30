function testLSC_DGS(n)


[uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, uI, vI, pI, infFlow] = dataOseen(n);

error = 1; step = 0;

while (error > 1e-6)
    [uh, vh, ph] = LSC_DGS(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig);
    
    [t1, t2, t3] = getResidualOseen(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig);
    
    error = sqrt(((norm(t1(:)))^2 + (norm(t2(:)))^2 + (norm(t3(:)))^2) / ((norm(f1h(:)))^2 + (norm(f2h(:)))^2 + (norm(gh(:)))^2))   
    step = step + 1;
end
% ph = ph(:); ph = ph - mean(ph); ph = reshape(ph, n, n);

[errVelL2, errVelH1, errVelInfi, errPreL2] = Error(uh, vh, ph, uI, vI, pI);
[errVelL2, errVelH1, errVelInfi, errPreL2, step]


uI = reshape(uI, n, n+1);  vI = reshape(vI, n+1, n); pI = reshape(pI, n, n);

% figure, surf(uh - uI); shading interp; colorbar;
% figure, surf(vh - vI); shading interp; colorbar;
% figure, surf(ph - pI); shading interp; colorbar;

end