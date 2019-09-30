function [r, r1, r2, r3] = assembleStiffnessMatrixInterfaceTermP0P0RT0_complex(NE, NGamma, psi, eitaGamma, brother, fracLength)

co1 = 4*ones(NGamma,1)./((2*psi-1)*eitaGamma); 
co2 = ones(NGamma,1)./eitaGamma;

r  = sparse(NE, NE);      r1 = sparse(NE, NGamma); 
r2 = sparse(NGamma, NE);  r3 = sparse(NGamma, NGamma);

for i = 1:NGamma
    j1 = brother(i,1); j2 = brother(i,3);
    
    h = fracLength(i);
    
    % if j1 is test 
    r(j1,j1) = r(j1,j1) + (co1(i)/4 + co2(i))*h;
    r(j1,j2) = r(j1,j2) + (co1(i)/4 - co2(i))*h;
    
    % if j2 is test 
    r(j2,j1) = r(j2,j1) + (co1(i)/4 - co2(i))*h;
    r(j2,j2) = r(j2,j2) + (co1(i)/4 + co2(i))*h;
    
    r1(j1,i) = r1(j1,i) - co1(i)/2*h;
    r1(j2,i) = r1(j2,i) - co1(i)/2*h;
    
    r2(i,j1) = r2(i,j1) - co1(i)/2*h;
    r2(i,j2) = r2(i,j2) - co1(i)/2*h;
    
    r3(i,i) = r3(i,i) + co1(i)*h;
end
end