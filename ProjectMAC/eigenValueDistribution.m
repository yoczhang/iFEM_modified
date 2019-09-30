%% 
% generate 4 largest and smallest eigenvalues of matrix 
% A = [-nu\Lap -\omega; \omega -nu\Lap];

k = 1;  
nu = [1e1 1 1e-1 1e-2 1e-3];
for i = 6:k+6
    for j = 1:length(nu)
        n = 2^i;
        [lamdaLarge, lamdaSmall] = test(n, nu(j));
        X = ['Viscosity is ', num2str(nu(j)), ', mesh size is ', num2str(n)];
        disp(X);
        disp('Large eigenvalues are '); disp(lamdaLarge);
        disp('Small eigenvalues are '); disp(lamdaSmall);
    end
end
