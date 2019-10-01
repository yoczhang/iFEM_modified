%% Vcycle Multigrid with Overweighting
% format short e
k = 1; infFlow = 1; nu = 1e-0; gamma = 1e-0;
ite = zeros(k+1,1); time = zeros(k+1,1); rate = zeros(k+1,1);
errVel = zeros(k+1,1); errPre = zeros(k+1,1); residual = cell(k+1,1);
for i = 6:k+6
    n = 2^i; level = i-1;
    [data, uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, uI, vI, pI] = dataSaddlePoint(n, nu, gamma);
    
    [uh, vh, ph, ite(i-5), time(i-5), residual{i-5}] = Vcycle(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, level, data.flowa, data.flowb, infFlow, nu, 'accurateSolution');
    ph = ph - mean(ph(:));
    
    ue = uh(:) - uI; ve = vh(:) - vI; pe = ph(:) - pI;
    uL2 = 1/n*norm(ue);  vL2 = 1/n*norm(ve);
    errVel(i-5) = sqrt(uL2^2 + vL2^2); errPre(i-5) = 1/n*norm(pe);
end
for i = 1:k+1
    for j = 4:length(residual{i})
        rate(i) = rate(i) + exp(log(residual{i}(j)/residual{i}(4)) / (j-3));
    end
    rate(i) = rate(i) / length(residual{i});
end
size = 2.^(6:k+6); size = size';
display('Table 1: Wcycle');
colname = {'#1/h', 'errVel', 'errPre', 'Ite', 'Rate', 'Time'};
disptable(colname,size,[], errVel,[], errPre,[], ite,[], rate,[], time,[]);
clear all


