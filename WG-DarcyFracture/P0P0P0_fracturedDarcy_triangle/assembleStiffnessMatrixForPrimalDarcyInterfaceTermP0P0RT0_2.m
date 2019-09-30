function [r1, r2, r3, r4, r5, r6, r7, r8] = assembleStiffnessMatrixForPrimalDarcyInterfaceTermP0P0RT0_2(NE1, ...
                                                            NE2, dof1D, brother, n, psi, eitaGamma, ...
                                                            Tb, basis_type)

r1 = sparse(NE1,NE1);   r2 = sparse(NE1,NE2);   r3 = sparse(NE2,NE1);
r4 = sparse(NE2,NE2);   r5 = sparse(dof1D,NE1); r6 = sparse(dof1D,NE2);
r7 = sparse(NE1,dof1D); r8 = sparse(NE2,dof1D);

co1 = 4/((2*psi-1)*eitaGamma); co2 = 1/eitaGamma;
h = 2 / n;
for k = 1:n
    j1 = brother(k,1); j2 = brother(k,2);
    r1(j1,j1) = r1(j1,j1) + (co1/4 + co2)*h;
    r2(j1,j2) = r2(j1,j2) + (co1/4 - co2)*h;
    r3(j2,j1) = r3(j2,j1) + (co1/4 - co2)*h;
    r4(j2,j2) = r4(j2,j2) + (co1/4 + co2)*h;
    
    beginPoint = -1+(k-1)*h; endPoint = -1+k*h; middle = (beginPoint + endPoint) / 2;
%     temp = - co1/2*h*localBasis1D(middle, 0, basis_type, beginPoint, endPoint);
    temp = - 0.5*co1*h*[0.5, 0.5];
    for beta = 1:basis_type+1
        r5(Tb(beta,k),j1) = r5(Tb(beta,k),j1) + temp(1,beta);
        r6(Tb(beta,k),j2) = r6(Tb(beta,k),j2) + temp(1,beta);
        
        r7(j1,Tb(beta,k)) = r7(j1,Tb(beta,k)) + temp(1,beta);
        r8(j2,Tb(beta,k)) = r8(j2,Tb(beta,k)) + temp(1,beta);
    end
end

end