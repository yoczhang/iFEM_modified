%%
% no numerical viscosity, take real viscosity value
fprintf('Real viscosity \n');
n = 128; viscosity = 1e-6; infinityNorm = 0; 
nu = viscosity;

directSolverForOseen(n, viscosity, infinityNorm, nu);

clear all

%% 
% take nu be numerical viscosity
fprintf('Numerical viscosity \n');
n = 128; viscosity = 1e-6; infinityNorm = 2;
nu = infinityNorm/(2*n);

directSolverForOseen(n, viscosity, infinityNorm, nu);

clear all