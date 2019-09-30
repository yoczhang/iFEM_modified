function [A, B] = modifyMatrix(Nu, Np, bdFlag, elem2dof, edge, A, B)

% First, find fixed points n
% Second, let A(n,:) = 0, A(:,n) = 0, A(n,n) = 1
N = Np;

isFixedDof = false(Nu,1);     
elem2edge = elem2dof(:,4:6)-N;
isDirichlet(elem2edge(bdFlag(:)==1)) = true;
isFixedDof(edge(isDirichlet,:)) = true;   % nodes of all D-edges
isFixedDof(N + find(isDirichlet')) = true;% dof on D-edges
fixedDof = find(isFixedDof);
ufreeDof = find(~isFixedDof); 

freeNode = [fixedDof; fixedDof+Nu];
bigA = [A B'; B sparse(N,N)];

bdidx = zeros(2*Nu+Np,1); 
bdidx(freeNode) = 1;
Tbd = spdiags(bdidx,0,2*Nu+Np,2*Nu+Np);
T = spdiags(1-bdidx,0,2*Nu+Np,2*Nu+Np);
bigAD = T*bigA*T + Tbd;

A = bigAD(1:2*Nu,1:2*Nu);
B = bigAD(2*Nu+1:end,1:2*Nu);

clear Tbd T bigAD

end