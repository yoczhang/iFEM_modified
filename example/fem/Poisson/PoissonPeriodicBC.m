%% RATE OF CONVERGENCE OF LINEAR ELEMENT FOR POISSON EQUATION
%
% This example is to show the rate of convergence of linear finite element
% approximation of the Poisson equation on the unit square with the
% following boundary conditions:
%
% - Non-empty Dirichlet boundary condition.
% - Pure Neumann boundary condition.
% - Robin boundary condition.
%
% The basis, data structure and numerical test is summarized in <a
% href="matlab:ifem Poissonfemrate">Poissonfemrate</a>.
%
% See also PoissonP2femrate, Poissonafemrate
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

clear variables
%% Setting
% 这里只验证很特殊的几个自由度的周期性边界条件, 不要更改网格大小和自由度排列
[node,elem] = squaremesh([0,1,0,1],0.25); 
mesh = struct('node',node,'elem',elem);
option.L0 = 0;
option.maxIt = 4;
option.printlevel = 1;
option.elemType = 'P1_periodic';

%% Non-empty Dirichlet boundary condition.
option.plotflag = 1;
pde = sincosdata;
mesh.bdFlag = setboundary(node,elem,'Dirichlet','(y==0) | (y==1)');
femPoisson(mesh,pde,option);

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% (2nd order) is observed. The 2nd order convergent rate between two
% discrete functions ||DuI-Duh|| is known as superconvergence.
%
% MGCG converges uniformly in all cases.
