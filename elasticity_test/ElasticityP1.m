function [solu, A, Mf]= ElasticityP1(node,elem,bdFlag,pde,option)
%% ELASTICITY  Conforming P1 elements discretization of linear elasticity equation
%
%   u = elasticity(node,elem,pde,bdEdge) use linear element to
%   approximate the displament u.
%
%       u = [ u1, u2]
%       -mu \Delta u - (lambda + mu)*grad(div(u)) = f in \Omega
%       Dirichlet boundary condition in Cartesian coordinates u = g_D, on \Gamma_D.
%       Dirichlet boundary condition in normal/tangential way u = g_S, on \Gamma_S.
%       Neumann boundary condition  sigma*n = g_N, on \Gamma_N.
%       Robin boundary condition    g_R*sigma + grad(\sigma)*n = g_N, on \Gamma_R.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('option','var'), option = []; end


%% Construct Data Structure
% important constants
N = size(node,1); NT = size(elem,1); Ndof = N;

%% Compute geometric quantities and gradient of local basis
[Dlambda,area] = gradbasis(node,elem);

%% Assemble stiffness matrix for (\nabla phi_i, \nabla phi_j)
A = sparse(N, N);

for i = 1:3
    for j = i:3
        Aij = pde.mu*(Dlambda(:,1,i).*Dlambda(:,1,j) + Dlambda(:,2,i).*Dlambda(:,2,j)).*area; 
        if j == i
            A = A + sparse(elem(:,i),elem(:,j),Aij, N, N);
        else
            A = A + sparse([elem(:,i);elem(:,j)], [elem(:,j);elem(:,i)],[Aij;Aij],N,N);
        end
    end
end
A = [A, sparse(N,N); sparse(N,N),A];

%% Assemble the matrix for ((lambda + mu)*div(phi_i), div(phi_j) )
B = sparse(2*N,2*N);
for i = 1:3
    for j = 1:3
        Bij11 = (pde.mu + pde.lambda)*Dlambda(:,1,i).*Dlambda(:,1,j).*area;
        Bij12 = (pde.mu + pde.lambda)*Dlambda(:,1,i).*Dlambda(:,2,j).*area;
        Bij21 = (pde.mu + pde.lambda)*Dlambda(:,2,i).*Dlambda(:,1,j).*area;
        Bij22 = (pde.mu + pde.lambda)*Dlambda(:,2,i).*Dlambda(:,2,j).*area;
        B = B + sparse(elem(:,i),   elem(:,j),     Bij11, 2*N, 2*N)...
              + sparse(elem(:,i),   elem(:,j) + N, Bij12, 2*N, 2*N)...
              + sparse(elem(:,i)+N, elem(:,j),     Bij21, 2*N, 2*N)...
              + sparse(elem(:,i)+N, elem(:,j) + N, Bij22, 2*N, 2*N);
    end
end
M = A + B;

%% Assemble right hand side by 4-points quadrature rule
f1 = zeros(N,1);
f2 = zeros(N,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 4;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end

if ~isempty(pde.f) 
    % quadrature points in the barycentric coordinate
    [lambda,weight] = quadpts(option.fquadorder); 
    phi = lambda;
    nQuad = length(weight);
    ft1 = zeros(NT, 3);
    ft2 = zeros(NT, 3);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        % function values at quadrature points
        fp = pde.f(pxy);
        % evaluate fp outside. 
        for j = 1:3
           ft1(:,j) = ft1(:,j) + fp(:,1).*phi(p,j)*weight(p);
           ft2(:,j) = ft2(:,j) + fp(:,2).*phi(p,j)*weight(p);
        end
    end
    ft1 = ft1.*repmat(area,1,3);
    ft2 = ft2.*repmat(area,1,3);
    f1 = accumarray(elem(:), ft1(:),[N,1]);
    f2 = accumarray(elem(:), ft2(:),[N,1]);
end
clear pxy ft1 ft2


%% Boundary condition

%% Initial check
if ~isfield(pde,'g_D'), pde.g_D = []; end
if ~isfield(pde,'g_N'), pde.g_N = []; end

% Find Dirichlet boundary nodes: fixedNode
fixedNode = []; freeNode = [];
if ~isempty(bdFlag) % find boundary edges and boundary nodes
    [fixedNode,bdEdge,isBdNode] = findboundary(elem,bdFlag);
    freeNode = find(~isBdNode);
end


%% Part 2: Find boundary edges and modify the right hand side b
% Find boundary edges: Neumann
Neumann = []; 
if ~isempty(bdFlag)  % bdFlag specifies different bd conditions
    Neumann = bdEdge;        
end

%% Neumann boundary condition
if  isnumeric(pde.g_N) && all(pde.g_N == 0)
    pde.g_N = [];
end
if ~isempty(Neumann) && ~isempty(pde.g_N)
    el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
    if ~isfield(option,'gNquadorder')
        option.gNquadorder = 2;   % default order exact for linear gN
    end
    [lambdagN,weightgN] = quadpts1(option.gNquadorder);
    phigN = lambdagN;                 % linear bases
    nQuadgN = size(lambdagN,1);
    ge1 = zeros(size(Neumann,1),2);
    ge2 = zeros(size(Neumann,1),2);
    for pp = 1:nQuadgN
        % quadrature points in the x-y coordinate
        ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
             + lambdagN(pp,2)*node(Neumann(:,2),:);
        gNp = pde.g_N(ppxy);
        for igN = 1:2
            ge1(:,igN) = ge1(:,igN) + weightgN(pp)*phigN(pp,igN)*gNp(:,1);
            ge2(:,igN) = ge2(:,igN) + weightgN(pp)*phigN(pp,igN)*gNp(:,2);
        end
    end
    ge1 = ge1.*repmat(el,1,2);
    ge2 = ge2.*repmat(el,1,2);  
    f1 = f1 + accumarray(Neumann(:), ge1(:),[Ndof,1]); 
    f2 = f2 + accumarray(Neumann(:), ge2(:),[Ndof,1]); 
end
% The case with non-empty Neumann edges but g_N=0 or g_N=[] corresponds to
% the zero flux boundary condition on Neumann edges and no modification of
% A,u,b is needed.       

% Dirichlet boundary condition
if isnumeric(pde.g_D) && all(pde.g_D == 0)   % zero g_D
    pde.g_D = [];
end

F = [f1;f2];
u = zeros(2*N,1);
ut = pde.g_D(node(fixedNode,:));
u(fixedNode) = ut(:,1); u(fixedNode+N) = ut(:,2);

freeDof = [freeNode;freeNode+N];

F = F - M*u;
u(freeDof) = M(freeDof, freeDof)\F(freeDof);

Mf = M(freeDof, freeDof);


%% Compute Du
u1 = u(1:N);
u2 = u(N+1:2*N);
du1dx =  u1(elem(:,1)).*Dlambda(:,1,1) + u1(elem(:,2)).*Dlambda(:,1,2) ...
       + u1(elem(:,3)).*Dlambda(:,1,3);
du1dy =  u1(elem(:,1)).*Dlambda(:,2,1) + u1(elem(:,2)).*Dlambda(:,2,2) ...
       + u1(elem(:,3)).*Dlambda(:,2,3);         
Du1 = [du1dx, du1dy];

du2dx =  u2(elem(:,1)).*Dlambda(:,1,1) + u2(elem(:,2)).*Dlambda(:,1,2) ...
      + u2(elem(:,3)).*Dlambda(:,1,3);
du2dy =  u2(elem(:,1)).*Dlambda(:,2,1) + u2(elem(:,2)).*Dlambda(:,2,2) ...
      + u2(elem(:,3)).*Dlambda(:,2,3);         
Du2 = [du2dx, du2dy];

%% Output
solu = struct('u',u,'u1',u1,'Du1',Du1,'u2',u2,'Du2',Du2);

end % elasticityP1
