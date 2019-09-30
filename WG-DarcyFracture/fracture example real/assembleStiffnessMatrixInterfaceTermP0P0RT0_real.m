function [r, r1, r2, r3] = assembleStiffnessMatrixInterfaceTermP0P0RT0_real(NE, NGamma, psi, eitaGamma, brother, fracLength)

co1 = 4/((2*psi-1)*eitaGamma); 
co2 = 1/eitaGamma;

r  = sparse(NE, NE);      r1 = sparse(NE, NGamma); 
r2 = sparse(NGamma, NE);  r3 = sparse(NGamma, NGamma);

for i = 1:NGamma
    j = brother(i,1); % fractured edge
    j1 = brother(i,2); j2 = brother(i,3); % matrix edge
    
    h = fracLength(i);
    
    % if j1 is test 
    r(j1,j1) = r(j1,j1) + (co1/4 + co2)*h;
    r(j1,j2) = r(j1,j2) + (co1/4 - co2)*h;
    
    % if j2 is test 
    r(j2,j1) = r(j2,j1) + (co1/4 - co2)*h;
    r(j2,j2) = r(j2,j2) + (co1/4 + co2)*h;
    
    r1(j1,j) = r1(j1,j) - co1/2*h;
    r1(j2,j) = r1(j2,j) - co1/2*h;
    
    r2(j,j1) = r2(j,j1) - co1/2*h;
    r2(j,j2) = r2(j,j2) - co1/2*h;
    
    r3(j,j) = r3(j,j) + co1*h;
end
end