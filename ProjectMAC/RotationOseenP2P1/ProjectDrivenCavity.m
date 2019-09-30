%% Setting
clear all; close all;
% [node,elem] = squaremesh([0,1,0,1],0.25);
node = [0,0; 1,0; 1,1; 0,1]; elem = [2,3,1; 4,1,3];
nu = 1e-2; pde = dataRotation2(nu); gamma = 1; 
bdFlag = setboundary(node,elem,'Dirichlet');

option.L0 = 6;
option.maxIt = 1;
option.printlevel = 1;
option.plotflag = 0;
option.rateflag = 0;

%% MG options
option.solver = 'mg';
option.smoothingstep = 2;
option.smootherbarSp = 'SGS';

%% P2-P1 element
display('P2-P1')
option.elemType = 'P2P1';
% option.smoothingStep = 3;
% option.smootherbarSp   = 'VCYCLE';
% option.smootherbarSpPara = 0.75;
[err,time,solver,eqn,u,p] = femStokes(node,elem,pde,bdFlag,option,nu,gamma);

for i = 1:option.L0
    [node,elem] = uniformrefine(node,elem);
end
for i = 1:option.maxIt-1
    [node,elem] = uniformrefine(node,elem);
end
[elem2dof,edge,bdDof] = dofP2(elem);
elem2 = [elem2dof(:,[1,6,5]); elem2dof(:,[2,4,6]); elem2dof(:,[3,5,4]); elem2dof(:,[4,5,6])];
node2 = [node; (node(edge(:,1),:)+node(edge(:,2),:))/2];

Nu = size(node,1) + size(edge,1); 
uh = u(1:Nu); vh = u(Nu+1:2*Nu); ph = p;
vel = sqrt(uh.^2 + vh.^2);
figure, trisurf(elem2, node2(:,1), node2(:,2), vel); shading interp; view(2); colorbar;
figure, trisurf(elem, node(:,1), node(:,2), ph); shading interp; view(2); colorbar;