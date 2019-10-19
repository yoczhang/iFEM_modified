%% Wcycle Multigrid with Overweighting
format short e
k = 3; infFlow = 1; nu = 1e-12;
ite = zeros(k+1,1); time = zeros(k+1,1); rate = zeros(k+1,1);
errVel = zeros(k+1,1); errPre = zeros(k+1,1); residual = cell(k+1,1);
for i = 6:k+6
    n = 2^i; level = i-1;
    [data, uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, uI, vI, pI] = dataOseen(n, nu);
    
    [uh, vh, ph, ite(i-5), time(i-5), residual{i-5}] = Wcycle(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, level, data.flowa, data.flowb, infFlow, nu, 'accurateSolution');
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


%% Defect Correction Procedure
format short e
k = 1; infFlow = 1; nu = 1e-12;
time = zeros(k+1,1); 
errVel = zeros(k+1,1); errPre = zeros(k+1,1);
for i = 6:k+6
    n = 2^i; level = i-1;
    [data, uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, uI, vI, pI] = dataOseen(n, nu);
    
    [uh, vh, ph, time(i-5)] = DefectionCorrection(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, level, data.flowa, data.flowb, infFlow, nu);
    ph = ph - mean(ph(:));
    
    ue = uh(:) - uI; ve = vh(:) - vI; pe = ph(:) - pI;
    uL2 = 1/n*norm(ue);  vL2 = 1/n*norm(ve); 
    errVel(i-5) = sqrt(uL2^2 + vL2^2); errPre(i-5) = 1/n*norm(pe);
end

size = 2.^(6:k+6); size = size';
display('Table 2: Defect Correction');
colname = {'#1/h', 'errVel', 'errPre', 'Time'};
disptable(colname,size,[], errVel,[], errPre,[], time,[]);
clear all



%% Driven cavity
format short e
k = 1; infFlow = 2; nu = 1e-6; rate = zeros(k+1,1);
ite = zeros(k+1,1); time = zeros(k+1,1); residual = cell(k+1,1);
for i = 6:k+6
    n = 2^i; level = i-1;
    [data, uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig] = dataOseenDrivenCavity(n);
    
    [uh, vh, ph, ite(i-5), time(i-5), residual{i-5}] = Wcycle(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, level, data.flowa, data.flowb, infFlow, nu, []);
end
 for i = 1:k+1
    for j = 4:length(residual{i})
        rate(i) = rate(i) + exp(log(residual{i}(j)/residual{i}(4)) / (j-3));
    end
    rate(i) = rate(i) / length(residual{i});
end
size = 2.^(6:k+6); size = size';
display('Table 2: Driven Cavity');
colname = {'#1/h', 'ite', 'Rate', 'Time'};
disptable(colname, size,[], ite,[], rate,[], time,[]);

% figure, surf(flipud(uh)); view(2); shading interp; colorbar;
% figure, surf(flipud(vh)); view(2); shading interp; colorbar;
% figure, surf(flipud(ph)); view(2); shading interp; colorbar;
clear all
