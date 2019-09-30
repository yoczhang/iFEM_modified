%%
% infinityNorm = ||u||_{infinity}
% assume that exact solution is zero, so infinityNorm==0
viscosity = 1e-12; infinityNorm = 1e-10;
size = [128; 256; 512];

[errVel, errPre] = firstStep(viscosity, infinityNorm, 'zeroInitial');
fprintf('Zero initial solution \n');
fprintf('#viscosity = %8.4e\n', viscosity);
colname = {'#1/h', 'errVel', 'errPre'};
disptable(colname, size,[], errVel,[], errPre,[]);


[errVel, errPre] = firstStep(viscosity, infinityNorm, 'randomInitial');
fprintf('Random initial solution \n');
fprintf('#viscosity = %8.4e\n', viscosity);
colname = {'#1/h', 'errVel', 'errPre'};
disptable(colname, size,[], errVel,[], errPre,[]);
