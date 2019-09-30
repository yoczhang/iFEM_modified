function [errVel, errPre] = firstStep(viscosity, infinityNorm, initialType)

errVel = zeros(3,1); errPre = zeros(3,1);
for i = 7:9
    n = 2^i; level = i-1; 
    if infinityNorm>0
        nu = infinityNorm/(2*n); 
    elseif infinityNorm==0
        nu = viscosity;
    end
%     exanType = 'accurateSolution';
    examType = [];
    [data, uh, vh, ph, f1h0, f2h0, gh, uLef, uRig, uTop, uBot, vLef, vRig, vTop, vBot, uI, vI, pI] = dataNavierStokes1(n, viscosity, initialType);
    [f1h, f2h] = modifyLoadVectorForLaplaceTerm(f1h0, f2h0, uLef, uRig, vTop, vBot, uTop, uBot, vLef, vRig, nu);
    [au, bu, av, bv] = formConvectionCoefficientMatrices(zeros(n,n+1), zeros(n+1,n), data.uexact, data.vexact, level);
    [f1h, f2h] = modifyLoadVectorForNonlinearTerm(f1h, f2h, uLef, uRig, uTop, uBot, vLef, vRig, vTop, vBot, au{level}, bu{level}, av{level}, bv{level});
    [uh, vh, ph] = WcycleNS(uh, vh, ph, f1h, f2h, gh, level, au, bu, av, bv, infinityNorm, viscosity, examType);

    ph = ph - mean(ph(:)); ue = uh(:) - uI; ve = vh(:) - vI; pe = ph(:) - pI;
    uL2 = 1/n*norm(ue);  vL2 = 1/n*norm(ve); 
    errVel(i-6) = sqrt(uL2^2 + vL2^2); errPre(i-6) = 1/n*norm(pe);
end

end