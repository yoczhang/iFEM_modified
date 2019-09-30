%% Setting
clear all; close all;
% [node,elem] = squaremesh([0,1,0,1],0.25);
node = [0,0; 1,0; 1,1; 0,1]; elem = [2,3,1; 4,1,3];
nu = 1e-3; pde = dataRotation(nu); gamma = 1;
bdFlag = setboundary(node,elem,'Dirichlet');

option.L0 = 0;
option.maxIt = 8;
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
femStokes(node,elem,pde,bdFlag,option,nu,gamma);