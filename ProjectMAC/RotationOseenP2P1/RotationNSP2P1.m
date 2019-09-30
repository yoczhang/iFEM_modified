function RotationNSP2P1(L0, nu)

%% Mesh
node = [0,0; 1,0; 1,1; 0,1]; elem = [2,3,1; 4,1,3];
bdFlag = setboundary(node,elem,'Dirichlet');
for k = 1:L0
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
end
[elem2dof, edge, bdDof] = dofP2(elem);
N = size(node,1);  NT = size(elem,1);  Nu = N+size(edge,1);  Np = N; 

data = dataRotation(nu);

%% Compute geometric quantities and gradient of local basis
[Dlambda,area] = gradbasis(node,elem);


%% Assemble stiffness matrix for Laplace operator
[A1, M] = formDiffusionMassMatrices(NT, elem2dof, Dlambda, area, node, elem, data.diffCo, data.omega, Nu);
A = [A1 -M; M A1]; clear A1 M


%% Assemble the matrix for divergence operator
B = divergenceOperator(Nu, Np, Dlambda, area, elem, elem2dof);


%% Assemble right hand side
[f1, f2] = formLoadVector(Nu, NT, elem2dof, area, node, elem, data.f);


%% Solve
isFixedDof = false(Nu,1);     
elem2edge = elem2dof(:,4:6)-N;
isDirichlet(elem2edge(bdFlag(:)==1)) = true;
isFixedDof(edge(isDirichlet,:)) = true;   % nodes of all D-edges
isFixedDof(N + find(isDirichlet')) = true;% dof on D-edges
fixedDof = find(isFixedDof);
ufreeDof = find(~isFixedDof);  

f = [f1; f2]; u1 = zeros(Nu,1); u2 = zeros(Nu,1); p = zeros(Np,1);

idx = (fixedDof > N);              % index of edge dof
uD = data.g_D(node(fixedDof(~idx),:));  % bd value at vertex dofs    
u1(fixedDof(~idx)) = uD(:,1); u2(fixedDof(~idx)) = uD(:,2);

bdEdgeIdx = fixedDof(idx)-N;
bdEdgeMid = (node(edge(bdEdgeIdx,1),:)+node(edge(bdEdgeIdx,2),:))/2;
uD = data.g_D(bdEdgeMid);         % bd values at middle points of edges
u1(fixedDof(idx)) = uD(:,1); u2(fixedDof(idx)) = uD(:,2);

u = [u1; u2]; % Dirichlet bd condition is built into u
f = f - A*u;  % bring affect of nonhomgenous Dirichlet bd condition to
g = zeros(Np,1); g = g - B*u;  % the right hand side
g = g - mean(g);         
f(fixedDof) = u1(fixedDof); f(fixedDof+Nu) = u2(fixedDof);
ufreeDof = [ufreeDof; ufreeDof+Nu]; pDof = (1:Np-1)';
bigA = [A, B'; ...
        B, sparse(Np,Np)];
bigF = [f; g];
bigu = [u; p];
bigFreeDof = [ufreeDof; 2*Nu+pDof];
bigu(bigFreeDof) = bigA(bigFreeDof,bigFreeDof)\bigF(bigFreeDof);
u = bigu(1:2*Nu);
p = bigu(2*Nu+1:end); 
if length(pDof) ~= Np % p is unique up to a constant
    % impose the condition int(p)=0
    c = sum(mean(p(elem),2).*area)/sum(area);
    p = p - c;
end


%% Figure

elem2 = [elem2dof(:,[1,6,5]); elem2dof(:,[2,4,6]); elem2dof(:,[3,5,4]); elem2dof(:,[4,5,6])];
node2 = [node; (node(edge(:,1),:)+node(edge(:,2),:))/2];
uI = data.exactu(node2); pI = data.exactp(node);
temp1 = norm(uI(:,1) - u(1:Nu)); temp2 = norm(uI(:,2) - u(Nu+1:2*Nu));
velL2 = sqrt(temp1^2 + temp2^2);
[size(node,1) velL2 norm(pI - p)]
% dofu = size(node2,1); 
% u1 = u(1:dofu,1); u2 = u(1+dofu:2*dofu,1); speed = sqrt(u1.^2 + u2.^2);
% 
% figure, trisurf(elem2, node2(:,1), node2(:,2), u(1:Nu)); shading interp;
% colorbar; xlabel('x'); ylabel('y'); view(2); axis equal; axis off;
% 
% figure, trisurf(elem2, node2(:,1), node2(:,2), u(Nu+1:2*Nu)); shading interp;
% colorbar; xlabel('x'); ylabel('y'); view(2); axis equal; axis off;
% 
% 
% figure, trisurf(elem, node(:,1), node(:,2), p); shading interp;
% colorbar; xlabel('x'); ylabel('y'); view(2); axis equal; axis off;
end