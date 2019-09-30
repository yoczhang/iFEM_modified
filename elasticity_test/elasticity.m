function [solu, A, Mf] = elasticity(node,elem,bdFlag,pde,option)
%% ELASTICITY  Conforming P1 elements discretization of linear elasticity equation
%
%   u = elasticity(node,elem,pde,bdEdge) use linear element to
%   approximate the displament u.
%
%       u = [ u1, u2]
%       -mu \Delta u - (lambda + mu)*grad(div(u)) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D.
%
% Created by Huayi Wei Monday, 27 June 2011.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


if ~exist('option','var'), option = []; end

N = size(node,1); NT = size(elem,1);

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
        
        Bij12_lambda = (pde.lambda)*Dlambda(:,1,i).*Dlambda(:,2,j).*area;
        Bij12_mu = (pde.mu)*Dlambda(:,2,i).*Dlambda(:,1,j).*area;
        
        Bij21_lambda = (pde.lambda)*Dlambda(:,2,i).*Dlambda(:,1,j).*area;
        Bij21_mu = (pde.mu)*Dlambda(:,1,i).*Dlambda(:,2,j).*area;
        
        Bij22 = (pde.mu + pde.lambda)*Dlambda(:,2,i).*Dlambda(:,2,j).*area;
        
        B = B + sparse(elem(:,i), elem(:,j), Bij11, 2*N, 2*N)...
            + sparse(elem(:,i), elem(:,j) + N, Bij12_lambda+Bij12_mu, 2*N, 2*N)...
            + sparse(elem(:,i)+N, elem(:,j), Bij21_lambda+Bij21_mu, 2*N, 2*N)...
            + sparse(elem(:,i)+N, elem(:,j) + N, Bij22, 2*N, 2*N);
    end
end
M = A + B;

%% Assemble right hand side by 4-points quadrature rule
f1 = zeros(N,1);
f2 = zeros(N,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 2;   % default order
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

%% Boundary condition
allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
Dirichlet = allEdge((bdFlag(:) == 1),:);
%---
% Dirichlet = Dirichlet(:);
% Neumann = allEdge((bdFlag(:) == 2) | (bdFlag(:) == 3),:);
% Neumann = Neumann(:);
% Dirichlet = setdiff(Dirichlet,intersect(Dirichlet,Neumann));
% close all
% showmesh(node,elem);
% findedge(node,Dirichlet,[],'noindex','LineWidth',2,'Color','r');
% findedge(node,Neumann,[],'noindex','LineWidth',2,'Color','k');
%---
isBdNode = false(N,1);
isBdNode(Dirichlet(:)) = true;
fixedNode = find(isBdNode);
% findnode(node,fixedNode(:),'noindex','color','r','MarkerSize',24)
% fixedNode = [fixedNode; fixedNode+N];
freeNode = find(~isBdNode);
% findnode(node,freeNode(:),'noindex','color','b','MarkerSize',24)
freeDof = [freeNode;freeNode+N];

u = zeros(2*N,1);
ut = pde.g_D(node(fixedNode,:));
u(fixedNode) = ut(:,1); u(fixedNode+N) = ut(:,2);

F = [f1;f2];
F = F - M(:,[fixedNode; fixedNode+N])*u([fixedNode; fixedNode+N]);

%% Solver
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