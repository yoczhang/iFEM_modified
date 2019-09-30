%%
% no numerical viscosity, take real viscosity value
fprintf('Real viscosity \n');
% viscosity = [1e-2 1e-4 1e-6 1e-8 1e-10]; 
viscosity = 1e-3;
infinityNorm = 0; n = length(viscosity);
errVel = zeros(3,n); errPre = zeros(3,n);
for k = 1:n
    for i = 6:8
        [errVel(i-5,k), errPre(i-5,k)] = directSolverForOseen(2^i, viscosity(k), infinityNorm);
    end
end

size = [64; 128; 256];
for k = 1:n
    fprintf('#viscosity = %8.4e\n', viscosity(k));
    colname = {'#1/h', 'errVel', 'errPre'};
    disptable(colname, size,[], errVel(:,k),[], errPre(:,k),[]);
end
clear all

%% 
% take nu be numerical viscosity
fprintf('Numerical viscosity \n');
viscosity = [1e-2 1e-4 1e-6 1e-8 1e-10]; 
infinityNorm = 1; n = length(viscosity);
errVel = zeros(3,n); errPre = zeros(3,n);
for k = 1:n
    for i = 6:8
        [errVel(i-5,k), errPre(i-5,k)] = directSolverForOseen(2^i, viscosity(k), infinityNorm);
    end
end

size = [64; 128; 256];
for k = 1:n
    fprintf('#viscosity = %8.4e\n', viscosity(k));
    colname = {'#1/h', 'errVel', 'errPre'};
    disptable(colname, size,[], errVel(:,k),[], errPre(:,k),[]);
end
clear all