function [A, B] = formDiffusionMassMatrices(NT, elem2dof, Dlambda, area, node, elem, diffCo, reacCo, Nu)

% generate sparse pattern
ii = zeros(36*NT,1); jj = zeros(36*NT,1); 
index = 0;
for i = 1:6
    for j = 1:6
        ii(index+1:index+NT) = double(elem2dof(:,i)); 
        jj(index+1:index+NT) = double(elem2dof(:,j));  
        index = index + NT;
    end
end
[lambda, w] = quadpts(9); nQuad = size(lambda,1);

phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
phi(:,4) = 4*lambda(:,2).*lambda(:,3);
phi(:,5) = 4*lambda(:,3).*lambda(:,1);
phi(:,6) = 4*lambda(:,1).*lambda(:,2);

% compute non-zeros
sA = zeros(36*NT,nQuad); sB = zeros(36*NT,nQuad);
for p = 1:nQuad
    % Dphi at quadrature points
    Dphip(:,:,6) = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
    Dphip(:,:,5) = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
    Dphip(:,:,4) = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
    Dphip(:,:,1) = (4*lambda(p,1)-1).*Dlambda(:,:,1);            
    Dphip(:,:,2) = (4*lambda(p,2)-1).*Dlambda(:,:,2);            
    Dphip(:,:,3) = (4*lambda(p,3)-1).*Dlambda(:,:,3);            
    index = 0;
    for i = 1:6
        for j = 1:6
            Aij = 0; Bij = 0;
            pxy = lambda(p,1)*node(elem(:,1),:) ...
                + lambda(p,2)*node(elem(:,2),:) ...
                + lambda(p,3)*node(elem(:,3),:);
            Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2).*diffCo(pxy);
            Bij = Bij + w(p)*phi(p,i).*phi(p,j).*reacCo(pxy);
            
            Aij = Aij.*area; Bij = Bij.*area;
            sA(index+1:index+NT,p) = Aij;
            sB(index+1:index+NT,p) = Bij;
            index = index + NT;
        end
    end
end
sA = sum(sA,2); sB = sum(sB,2);
% assemble the matrix
A = sparse(ii,jj,sA,Nu,Nu); B = sparse(ii,jj,sB,Nu,Nu);
clear Aij ii jj sA

end