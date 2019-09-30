function [u1, u2, f1, f2, g] = modifyDirichletBoundaryCondition(Nu, N, node, elem2dof, bdFlag, edge, A, B, f1, f2, g_D)

Np = N;
isFixedDof = false(Nu,1);     
elem2edge = elem2dof(:,4:6)-N;
isDirichlet(elem2edge(bdFlag(:)==1)) = true;
isFixedDof(edge(isDirichlet,:)) = true;   % nodes of all D-edges
isFixedDof(N + find(isDirichlet')) = true;% dof on D-edges
fixedDof = find(isFixedDof);
% ufreeDof = find(~isFixedDof);  

f = [f1; f2]; u1 = zeros(Nu,1); u2 = zeros(Nu,1); 

idx = (fixedDof > N);              % index of edge dof
uD = g_D(node(fixedDof(~idx),:));  % bd value at vertex dofs    
u1(fixedDof(~idx)) = uD(:,1); u2(fixedDof(~idx)) = uD(:,2);

bdEdgeIdx = fixedDof(idx)-N;
bdEdgeMid = (node(edge(bdEdgeIdx,1),:)+node(edge(bdEdgeIdx,2),:))/2;
uD =g_D(bdEdgeMid);         % bd values at middle points of edges
u1(fixedDof(idx)) = uD(:,1); u2(fixedDof(idx)) = uD(:,2);

u = [u1; u2]; % Dirichlet bd condition is built into u
f = f - A*u;  % bring affect of nonhomgenous Dirichlet bd condition to
g = zeros(Np,1); g = g - B*u;  % the right hand side
g = g - mean(g);         
f(fixedDof) = u1(fixedDof); f(fixedDof+Nu) = u2(fixedDof);
f1 = f(1:Nu); f2 = f(Nu+1:end);

end