%% MAC for SaddlePoint

close all
clear

k = 4; 
nn = 5;
mu = 1e-4; 
gamma = 1e-2;

errVelL2 = zeros(k+1,1); errVelH1 = zeros(k+1,1); errVelInfi = zeros(k+1,1);
errPreL2 = zeros(k+1,1);
ite = zeros(k+1,1); time = zeros(k+1,1);

tol = 1e-12;  maxit = 15; 
% errVel = zeros(k+1,1); errPre = zeros(k+1,1); residual = cell(k+1,1);
for i = nn:nn+k
    n = 2^i; level = i;
    %- generate the init values
    [uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, uI, vI, pI, h, width] = gmres_dataSaddlePoint(n, mu, gamma);
        %> actually, we only the f1h, f2h, gh
        
    bh = [f1h(:); f2h(:); gh(:)];
    Uh = gmres(@as_gmres_afun, bh, 10, tol, maxit, @as_gmres_mfun);
    
    
    %check_results(uh, vh, uI, vI)
    
    [errVelL2(i-nn+1), errVelH1(i-nn+1), errVelInfi(i-nn+1), errPreL2(i-nn+1)] = Error(uh, vh, ph, uI, vI, pI);
    
%     ue = uh(:) - uI; ve = vh(:) - vI; pe = ph(:) - pI;
%     uL2 = 1/n*norm(ue);  vL2 = 1/n*norm(ve);
%     errVel(i-nn+1) = sqrt(uL2^2 + vL2^2); errPre(i-nn+1) = 1/n*norm(pe);
end


size = 2.^(nn:nn+k); size = size';
disp('Table 1: Wcycle');
colname = {'#1/h', 'errVelL2', 'errVelH1', 'errVelInfi', 'errPreL2', 'Ite', 'Time'};
disptable(colname,size,[], errVelL2,[], errVelH1,[], errVelInfi,[], errPreL2,[], ite,[], time,[]);

showfigure = false;
if showfigure
    figure; showrate(width./size, errVelL2);
    h1 = legend('$ \frac {|| \mathbf{u}_I - \mathbf{u}_h ||} {|| \mathbf{u}_I ||}$','Location','southeast');
    set(h1,'Interpreter','latex')

    figure; showrate(width./size, errVelH1);
    h2 = legend('$|| \nabla(\mathbf{u}_I - \mathbf{u}_h) ||$','Location','southeast');
    set(h2,'Interpreter','latex')

    figure; showrate(width./size, errVelInfi);
    h3 = legend('$|| \mathbf{u}_I - \mathbf{u}_h ||_{\infty}$','Location','southeast');
    set(h3,'Interpreter','latex')

    figure; showrate(width./size, errPreL2);
    h4 = legend('$ \frac {|| p_I - p_h ||} {|| p_I||}$','Location','southeast');
    set(h4,'Interpreter','latex')
end



