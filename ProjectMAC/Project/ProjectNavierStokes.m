%% Analytical solution
viscosity = 1e-12; 
errVel = zeros(3,1); errPre = zeros(3,1);
for i = 6:8
    [errVel(i-5,1), errPre(i-5,1)] = subProjectNS(i, viscosity);
end

size = [128; 256; 512];
fprintf('#viscosity = %8.4e\n', viscosity);
colname = {'#1/h', 'errVel', 'errPre'};
disptable(colname, size,[], errVel,[], errPre,[]);

%% Driven cavity
% viscosity = 1e-6; infinityNorm = 2;
% [errVel, errPre] = subProjectNS(8, viscosity, infinityNorm);

