function [r1, r2, r3, r4, r5, r6, r7, r8] = assembleStiffnessMatrixForPrimalDarcyInterfaceTermP1P1RT1_2(NE1, ...
                                                            NE2, dof1D, brother, n, psi, eitaGamma, ...
                                                            Tb, basis_type, node1, edge1, node2, edge2, gcrl, gprl)

r1 = sparse(2*NE1,2*NE1); r2 = sparse(2*NE1,2*NE2); r3 = sparse(2*NE2,2*NE1);
r4 = sparse(2*NE2,2*NE2); r5 = sparse(dof1D,2*NE1); r6 = sparse(dof1D,2*NE2);
r7 = sparse(2*NE1,dof1D); r8 = sparse(2*NE2,dof1D);
h = 2 / n; % length of one edge
co1 = 4/((2*psi-1)*eitaGamma); co2 = 1/eitaGamma; dof = basis_type + 1;
for k = 1:n
    j1 = brother(k,1); j2 = brother(k,2); middle = -1+h/2+(k-1)*h;
    bp1 = node1(edge1(j1,1),:); ep1 = node1(edge1(j1,2),:); bp2 = node2(edge2(j2,1),:); ep2 = node2(edge2(j2,2),:);
%     [gcll,gpll] = generate_Gauss_local_1D(gcrl, gprl, [-1+(k-1)*h, 0], [-1+k*h, 0], 0);
    [gcll,gpll] = generate_Gauss_local_1D(gcrl, gprl, bp1, ep1, 0);
    
    temp1 = zeros(2,2);   temp2 = zeros(2,2);   temp3 = zeros(2,2);   temp4 = zeros(2,2);
    temp5 = zeros(2,dof); temp6 = zeros(2,dof); temp7 = zeros(dof,2); temp8 = zeros(dof ,2);
    for j = 1:length(gcll)
        temp1 = temp1 + gcll(j) * basis_1(gpll(j), middle, bp1, ep1)' * basis_1(gpll(j), middle, bp1, ep1);
        temp2 = temp2 + gcll(j) * basis_2(gpll(j), middle, bp2, ep2)' * basis_1(gpll(j), middle, bp1, ep1);
        temp3 = temp3 + gcll(j) * basis_1(gpll(j), middle, bp1, ep1)' * basis_2(gpll(j), middle, bp2, ep2);
        temp4 = temp4 + gcll(j) * basis_2(gpll(j), middle, bp2, ep2)' * basis_2(gpll(j), middle, bp2, ep2);
        
        temp5 = temp5 + gcll(j) * basis_1(gpll(j), middle, bp1, ep1)' * ...
                                  localBasis1D(gpll(j), 0, basis_type, [-1+(k-1)*h, 0], [-1+k*h, 0]);
        temp6 = temp6 + gcll(j) * basis_2(gpll(j), middle, bp2, ep2)' * ...
                                  localBasis1D(gpll(j), 0, basis_type, [-1+(k-1)*h, 0], [-1+k*h, 0]);
        temp7 = temp7 + gcll(j) * localBasis1D(gpll(j), 0, basis_type, [-1+(k-1)*h, 0], [-1+k*h, 0])' * ...
                                  basis_1(gpll(j), middle, bp1, ep1);
        temp8 = temp8 + gcll(j) * localBasis1D(gpll(j), 0, basis_type, [-1+(k-1)*h, 0], [-1+k*h, 0])' * ...
                                  basis_2(gpll(j), middle, bp2, ep2);                              
    end
    used1 = [2*j1-1,2*j1]; used2 = [2*j2-1,2*j2];
    for alpha = 1:2
        for beta = 1:2
            r1(used1(beta),used1(alpha)) = r1(used1(beta),used1(alpha)) + temp1(alpha,beta);
            r2(used1(beta),used2(alpha)) = r2(used1(beta),used2(alpha)) + temp2(alpha,beta);
            r3(used2(beta),used1(alpha)) = r3(used2(beta),used1(alpha)) + temp3(alpha,beta);
            r4(used2(beta),used2(alpha)) = r4(used2(beta),used2(alpha)) + temp4(alpha,beta);
        end
    end
    for alpha = 1:2
        for beta = 1:dof
            r5(Tb(beta,k),used1(alpha)) = r5(Tb(beta,k),used1(alpha)) + temp5(alpha,beta);
            r6(Tb(beta,k),used2(alpha)) = r6(Tb(beta,k),used2(alpha)) + temp6(alpha,beta);
        end
    end
    for alpha = 1:dof
        for beta = 1:2
            r7(used1(beta),Tb(alpha,k)) = r7(used1(beta),Tb(alpha,k)) + temp7(alpha,beta);
            r8(used2(beta),Tb(alpha,k)) = r8(used2(beta),Tb(alpha,k)) + temp8(alpha,beta);
        end
    end   
end
r1 = (co1/4+co2)*r1; r2 = (co1/4-co2)*r2; r3 = (co1/4-co2)*r3; r4 = (co1/4+co2)*r4;
r5 = -co1/2*r5;      r6 = -co1/2*r6;      r7 = -co1/2*r7;      r8 = -co1/2*r8;
    function s = basis_1(x, xc, begin_1, ending_1) % domain 1
        s = [1, (x-xc)/((ending_1(1)-begin_1(1))/2)];
    end
    function s = basis_2(x, xc, begin_2, ending_2) % domain 2
        s = [1, (x-xc)/((ending_2(1)-begin_2(1))/2)];
    end
end