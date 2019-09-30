function [r1, r2, r3, r4, r5, r6, r7, r8] = assembleStiffnessMatrixForPrimalDarcyInterfaceTermP0P0RT0P0DG_2(NE1, ...
                                                            NE2, brother, n, psi, eitaGamma)

r1 = sparse(NE1,NE1); r2 = sparse(NE1,NE2); r3 = sparse(NE2,NE1);
r4 = sparse(NE2,NE2); r5 = sparse(n,NE1); r6 = sparse(n,NE2);
r7 = sparse(NE1,n); r8 = sparse(NE2,n);

co1 = 4/((2*psi-1)*eitaGamma); co2 = 1/eitaGamma;
h = 1 / n;
for k = 1:n
    j1 = brother(k,1); j2 = brother(k,2);
    r1(j1,j1) = r1(j1,j1) + (co1/4 + co2)*h;
    r2(j1,j2) = r2(j1,j2) + (co1/4 - co2)*h;
    r3(j2,j1) = r3(j2,j1) + (co1/4 - co2)*h;
    r4(j2,j2) = r4(j2,j2) + (co1/4 + co2)*h;
    
        r5(k,j1) = r5(k,j1) - 0.5*co1*h;
        r6(k,j2) = r6(k,j2) - 0.5*co1*h;
        
        r7(j1,k) = r7(j1,k) - 0.5*co1*h;
        r8(j2,k) = r8(j2,k) - 0.5*co1*h;
end