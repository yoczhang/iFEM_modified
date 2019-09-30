function [errVel, errPre] = subProjectNS(i, viscosity)
%% Analytical solution 
n = 2^i; level = i-1; 
% X = ['Mesh size is ', num2str(n)]; disp(X);
% fprintf('The first part\n'); ite = 1; 

%% Stokes
% % Firstly, set initial value u = 0, that is, solve an Stokes problem
% % in this case, au, bu, av, bv are all zeros
% % and using real viscosity coefficient
% [data, uh, vh, ph, f1h0, f2h0, gh, uLef, uRig, uTop, uBot, vLef, vRig, vTop, vBot, uI, vI, pI] = dataNavierStokes(n, viscosity);
% [f1h, f2h] = modifyLoadVectorForLaplaceTerm(f1h0, f2h0, uLef, uRig, vTop, vBot, uTop, uBot, vLef, vRig, viscosity);
% [au, bu, av, bv] = formConvectionCoefficientMatrices(zeros(n,n+1), zeros(n+1,n), data.emptyFunc, data.emptyFunc, level);
% [uh, vh, ph] = WcycleNS(uh, vh, ph, f1h, f2h, gh, level, au, bu, av, bv, 0, viscosity, 'accurateSolution');
% uh3 = uh; vh3 = vh; Error = sqrt(norm(uh-uh3,'fro')^2 + norm(vh-vh3,'fro')^2);
% fprintf('#Outer step: %8.0u, OUTER err = %8.4e\n',ite,Error); clear uh3 vh3
% 
% 
% %% Oseen loop
% % we have obtain initial solutions from Stokes equation with real viscosity
% % now using w-cycle multigrid with overweighting 
% Error = 1; 
% % viscosity = 2*viscosity;
% % infinityNorm = max(abs(uh(:))) + max(abs(vh(:)));
% % [data, uh, vh, ph, f1h, f2h, gh, uLef, uRig, uTop, uBot, vLef, vRig, vTop, vBot, uI, vI, pI] = dataNavierStokes(n, viscosity);
% % [f1h, f2h] = modifyLoadVectorForLaplaceTerm(f1h, f2h, uLef, uRig, vTop, vBot, uTop, uBot, vLef, vRig, viscosity);
% while (Error>1e-6)
%     [au, bu, av, bv] = formConvectionCoefficientMatrices(uh, vh, data.uexact, data.vexact, level);
%     [f1h, f2h] = modifyLoadVectorForLaplaceTerm(f1h0, f2h0, uLef, uRig, vTop, vBot, uTop, uBot, vLef, vRig, viscosity);
%     [f1h, f2h] = modifyLoadVectorForNonlinearTerm(f1h, f2h, uLef, uRig, uTop, uBot, vLef, vRig, vTop, vBot, au{level}, bu{level}, av{level}, bv{level});
% %     infinityNorm = max(abs(uh(:))) + max(abs(vh(:)));
%     infinityNorm = 0;
%     [uh2, vh2, ph2] = WcycleNS(uh, vh, ph, f1h, f2h, gh, level, au, bu, av, bv, infinityNorm, viscosity, 'accurateSolution');
% 
%     Error = sqrt(norm(uh-uh2,'fro')^2 + norm(vh-vh2,'fro')^2);
%     ite = ite + 1;  uh = uh2; vh = vh2; ph = ph2; 
%     fprintf('#Outer step: %8.0u, OUTER err = %8.4e\n', ite, Error);  
% end
% 
% % Loop
% % Error = 1; viscosity = viscosity/2; 
% % [~, ~, ~, ~, f1h, f2h, gh, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = dataNavierStokes(n, viscosity);
% % while (Error>1e-6)
% %     [au, bu, av, bv] = formConvectionCoefficientMatrices(uh, vh, data.uexact, data.vexact, level);
% %     [f1h, f2h] = modifyLoadVectorForNonlinearTerm(f1h, f2h, uLef, uRig, uTop, uBot, vLef, vRig, vTop, vBot, au{level}, bu{level}, av{level}, bv{level});
% %     infinityNorm = max(abs(uh(:))) + max(abs(vh(:)));
% % %     infinityNorm = 0;
% %     [uh2, vh2, ph2] = WcycleNS(uh, vh, ph, f1h, f2h, gh, level, au, bu, av, bv, infinityNorm, viscosity, 'accurateSolution');
% % 
% %     Error = sqrt(norm(uh-uh2,'fro')^2 + norm(vh-vh2,'fro')^2);
% %     ite = ite + 1;  uh = uh2; vh = vh2; ph = ph2; 
% %     fprintf('#Outer step: %8.0u, OUTER err = %8.4e\n', ite, Error);  
% % end
% 
% ph = ph - mean(ph(:)); ue = uh(:) - uI; ve = vh(:) - vI; pe = ph(:) - pI;
% uL2 = 1/n*norm(ue);  vL2 = 1/n*norm(ve); 
% errVel = sqrt(uL2^2 + vL2^2); errPre = 1/n*norm(pe);
% 
% X = ['The viscosity is ', num2str(viscosity)]; disp(X);
% X = ['The errors of velocity and pressure are ', num2str([errVel, errPre])]; disp(X);
% clear all



%% Oseen with given convection coefficient
% refer to paper: A multigrid solver based on... Long Chen, Xiaozhe Hu...

[data, uh, vh, ph, f1h, f2h, gh, uLef, uRig, uTop, uBot, vLef, vRig, vTop, vBot, uI, vI, pI] = dataNavierStokes(n, viscosity);
[au, bu, av, bv] = formConvectionCoefficientMatrices(zeros(n,n+1), zeros(n+1,n), data.flowa, data.flowb, level);
infinityNorm = max([abs(au{level}(:)); abs(bu{level}(:)); abs(av{level}(:)); abs(bv{level}(:))]);
h = 1/n; nu = infinityNorm*h/2;
[f1h, f2h] = modifyLoadVectorForLaplaceTerm(f1h, f2h, uLef, uRig, vTop, vBot, uTop, uBot, vLef, vRig, nu);
[f1h, f2h] = modifyLoadVectorForNonlinearTerm(f1h, f2h, uLef, uRig, uTop, uBot, vLef, vRig, vTop, vBot, au{level}, bu{level}, av{level}, bv{level});
[uh, vh, ph] = WcycleNS(uh, vh, ph, f1h, f2h, gh, level, au, bu, av, bv, infinityNorm, viscosity, 'accurateSolution');

ue = uh(:) - uI; ve = vh(:) - vI; pe = ph(:) - pI;
uL2 = 1/n*norm(ue);  vL2 = 1/n*norm(ve); 
errVel = sqrt(uL2^2 + vL2^2); errPre = 1/n*norm(pe);

% figure, surf(flipud(uh)); view(2); shading interp; colorbar;
% figure, surf(flipud(vh)); view(2); shading interp; colorbar;
% figure, surf(flipud(ph)); view(2); shading interp; colorbar;
% X = ['The viscosity is ', num2str(viscosity)]; disp(X);
% X = ['The errors of velocity and pressure are ', num2str([errVel, errPre])]; disp(X);
end