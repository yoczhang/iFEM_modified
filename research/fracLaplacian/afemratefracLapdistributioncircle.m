close all; clear all;

%% Setting of the problem
global s
pde = fracLapdata6; 
option.theta = 0.3;
option.estType = 'star';
option.maxIt = 18;
option.maxN = 5e4;
option.solver = 'mg';
option.tol = 1e-6;
[node,elem] = squaremesh([-1,1,-1,1],0.5);
bdFlag = setboundary(node,elem,'Dirichlet');

% %% s = 0.2
% s = 0.2;
% afemfracLap(node,elem,pde,bdFlag,option);
% 
% %% s = 0.4
% s = 0.4;
% afemfracLap(node,elem,pde,bdFlag,option);

%% s = 0.6
s = 0.6;
afemfracLap(node,elem,pde,bdFlag,option);

%% s = 0.8
s = 0.8;
afemfracLap(node,elem,pde,bdFlag,option);