function [I1,I2,M,D,rhs] = ...
    assembleInterfaceMatrix_FVMmatrix(node,edge,lGamma,eipsilon,psi,eitaGamma,...
    NE,inteNE,edgeLength, ...
    interfacesEdges_domainOne,interfacesEdges_domainTwo,pGamma,fGamma)

co1 = 4/((2*psi-1)*eitaGamma); co2 = 1/eitaGamma;
I1 = sparse(NE,NE); I2 = sparse(NE,inteNE); M = sparse(inteNE,inteNE); rhs = zeros(inteNE,1);
alpha = 2*lGamma*eipsilon./edgeLength(interfacesEdges_domainOne);
t = [0, sqrt(21)/7, -sqrt(21)/7, 1, -1];  weight = [32/45, 49/90, 49/90, 1/10, 1/10];

for i = 1:inteNE
    j1 = interfacesEdges_domainOne(i); 
    j2 = interfacesEdges_domainTwo(i); 
    I1(j1,j1) = I1(j1,j1) + (co1/4 + co2) * edgeLength(j1); 
    I1(j2,j1) = I1(j2,j1) + (co1/4 - co2) * edgeLength(j1);
    I1(j1,j2) = I1(j1,j2) + (co1/4 - co2) * edgeLength(j1); 
    I1(j2,j2) = I1(j2,j2) + (co1/4 + co2) * edgeLength(j1);
    
    I2(j1,i) = I2(j1,i) - co1/2  * edgeLength(j1); 
    I2(j2,i) = I2(j2,i) - co1/2  * edgeLength(j1);
    
    M(i,i) = M(i,i) + co1 * edgeLength(j1);
    for j = 1:length(t)
        left = node(edge(j1,1),:); righ = node(edge(j1,2),:); pxy = (righ - left)/2 * t(j) + (righ + left)/2;
        rhs(i,1) = rhs(i,1) + lGamma * edgeLength(j1)/2 * weight(j) * fGamma(pxy(1));
    end
end

temp = alpha(2:end-1).*alpha(3:end) ./ (alpha(2:end-1) + alpha(3:end)) + alpha(2:end-1).*alpha(1:end-2) ./ (alpha(2:end-1) + alpha(1:end-2));
rhs(1) = rhs(1) + alpha(1)*pGamma(0); rhs(end) = rhs(end) + alpha(end)*pGamma(1);
ii = [alpha(1) + alpha(1)*alpha(2) / (alpha(1)+alpha(2)); temp; alpha(end) + alpha(end-1)*alpha(end)/(alpha(end-1)+alpha(end))];
D = spdiags(ii,0,inteNE,inteNE);
for i = 1:inteNE-1
    D(i,i+1) = D(i,i+1) - alpha(i)*alpha(i+1)/(alpha(i)+alpha(i+1));
    D(i+1,i) = D(i+1,i) - alpha(i)*alpha(i+1)/(alpha(i)+alpha(i+1));
end

end