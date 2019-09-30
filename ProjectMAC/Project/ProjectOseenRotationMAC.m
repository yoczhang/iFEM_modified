%%
format short e
k = 2; nu = 1e-0;
ite = zeros(k+1,1); time = zeros(k+1,1); rate = zeros(k+1,1);
errVelL2 = zeros(k+1,1); errVelH1 = zeros(k+1,1);
errVelInfi = zeros(k+1,1); errPreL2 = zeros(k+1,1); residual = cell(k+1,1);
for i = 6:k+6
    n = 2^i; level = i-1;
    [data, uh, vh, ph, f1h, f2h, gh, uI, vI, pI, width, nu] = dataOseenRotation(n, nu);
    
    [uh, vh, ph, ite(i-5), time(i-5), residual{i-5}] = VcycleRotation(uh, vh, ph, f1h, f2h, gh, data.omega, level, width, nu);
    
    [errVelL2(i-5), errVelH1(i-5), errVelInfi(i-5), errPreL2(i-5)] = Error(uh, vh, ph, uI, vI, pI);
end

for i = 1:k+1
    for j = 4:length(residual{i})
        rate(i) = rate(i) + exp(log(residual{i}(j)/residual{i}(4)) / (j-3));
    end
    rate(i) = rate(i) / length(residual{i});
end
size = 2.^(6:k+6); size = size';
display('Table 1: Vcycle');
colname = {       '#1/h', 'errVelL2',  'errVelH1',  'errVelInfi',  'errPreL2',  'Ite',  'Rate',  'Time'};
disptable(colname,size,[], errVelL2,[], errVelH1,[], errVelInfi,[], errPreL2,[], ite,[], rate,[], time,[]);
clear all