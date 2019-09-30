function [r, r1, r2, r3, NGamma] = assembleStiffnessMatrixInterfaceTermP0P0RT0_geiger(NE, h, frac1, frac2, frac3, frac4, frac5, frac6, psi, eitaGamma)

r = sparse(NE, NE);
frac = [frac1; frac2; frac3; frac4; frac5; frac6];
co1 = 4/((2*psi-1)*eitaGamma); co2 = 1/eitaGamma;

NGamma = size(frac,1);
r1 = sparse(NE, NGamma); r2 = sparse(NGamma, NE);
r3 = sparse(NGamma, NGamma);

for i = 1:NGamma
    j1 = frac(i,1); j2 = frac(i,3);
    
    % if j1 is test 
    r(j1,j1) = r(j1,j1) + (co1/4 + co2)*h;
    r(j1,j2) = r(j1,j2) + (co1/4 - co2)*h;
    
    % if j2 is test 
    r(j2,j1) = r(j2,j1) + (co1/4 - co2)*h;
    r(j2,j2) = r(j2,j2) + (co1/4 + co2)*h;
    
    r1(j1,i) = r1(j1,i) - co1/2*h;
    r1(j2,i) = r1(j2,i) - co1/2*h;
    
    r2(i,j1) = r2(i,j1) - co1/2*h;
    r2(i,j2) = r2(i,j2) - co1/2*h;
    
    r3(i,i) = r3(i,i) + co1*h;
end



end