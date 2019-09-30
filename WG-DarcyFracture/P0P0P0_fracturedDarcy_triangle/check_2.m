%%
n = 6; eipsilon = 1;
error0 = zeros(n,1); errorb = zeros(n,1); H1p = zeros(n,1);
H1pGamma = zeros(n,1); maxpGamma = zeros(n,1); dof = zeros(n,1);
for i = 1:n
    [error0(i), errorb(i), H1p(i), H1pGamma(i), maxpGamma(i), dof(i)] = FractureDarcyP0P0RT0_P0DG_2(eipsilon, i);
end
error0Rate = zeros(n,1); errorbRate = zeros(n,1); H1pRate = zeros(n,1);
H1pGammaRate = zeros(n,1); maxpGammaRate = zeros(n,1); digits(3);
for i = 2:n
    error0Rate(i) = vpa( log(error0(i-1)/error0(i))/log(2));
    errorbRate(i) = vpa( log(errorb(i-1)/errorb(i))/log(2));
    H1pRate(i) = vpa( log(H1p(i-1)/H1p(i))/log(2));
    H1pGammaRate(i) = vpa( log(H1pGamma(i-1)/H1pGamma(i))/log(2));
    maxpGammaRate(i) = vpa( log(maxpGamma(i-1)/maxpGamma(i))/log(2));
end
display('Table 1: error when eipsilon is 1');
colname = {'#DOF',         'error0',   'errorb', 'H1p',  'H1pGamma', 'maxpGamma'};
disptable(colname,dof(3:n),[], error0(3:n),'%0.4e', errorb(3:n),'%0.4e', H1p(3:n), '%0.4e',  H1pGamma(3:n),'%0.4e', maxpGamma(3:n),'%0.4e');
rate = [sum(error0Rate(3:n))/(n-2), sum(errorbRate(3:n))/(n-2), sum(H1pRate(3:n))/(n-2), sum(H1pGammaRate(3:n))/(n-2), sum(maxpGammaRate(3:n))/(n-2)];
for i = 1:5, if i<5, fprintf(' %0.2f  ', rate(i)); else fprintf(' %0.2f  \n', rate(i)); end; end
clear all

%%
n = 6; eipsilon = 1000;
error0 = zeros(n,1); errorb = zeros(n,1); H1p = zeros(n,1);
H1pGamma = zeros(n,1); maxpGamma = zeros(n,1); dof = zeros(n,1);
for i = 1:n
    [error0(i), errorb(i), H1p(i), H1pGamma(i), maxpGamma(i), dof(i)] = FractureDarcyP0P0RT0_P0DG_2(eipsilon, i);
end
error0Rate = zeros(n,1); errorbRate = zeros(n,1); H1pRate = zeros(n,1);
H1pGammaRate = zeros(n,1); maxpGammaRate = zeros(n,1); digits(3);
for i = 2:n
    error0Rate(i) = vpa( log(error0(i-1)/error0(i))/log(2));
    errorbRate(i) = vpa( log(errorb(i-1)/errorb(i))/log(2));
    H1pRate(i) = vpa( log(H1p(i-1)/H1p(i))/log(2));
    H1pGammaRate(i) = vpa( log(H1pGamma(i-1)/H1pGamma(i))/log(2));
    maxpGammaRate(i) = vpa( log(maxpGamma(i-1)/maxpGamma(i))/log(2));
end
display('Table 2: error when eipsilon is 1000');
colname = {'#DOF',         'error0',   'errorb', 'H1p',  'H1pGamma', 'maxpGamma'};
disptable(colname,dof(3:n),[], error0(3:n),'%0.4e', errorb(3:n),'%0.4e', H1p(3:n), '%0.4e',  H1pGamma(3:n),'%0.4e', maxpGamma(3:n),'%0.4e');
rate = [sum(error0Rate(3:n))/(n-2), sum(errorbRate(3:n))/(n-2), sum(H1pRate(3:n))/(n-2), sum(H1pGammaRate(3:n))/(n-2), sum(maxpGammaRate(3:n))/(n-2)];
for i = 1:5, if i<5, fprintf(' %0.2f  ', rate(i)); else fprintf(' %0.2f  \n', rate(i)); end; end
clear all

%%
n = 6; eipsilon = 0.001;
error0 = zeros(n,1); errorb = zeros(n,1); H1p = zeros(n,1);
H1pGamma = zeros(n,1); maxpGamma = zeros(n,1); dof = zeros(n,1);
for i = 1:n
    [error0(i), errorb(i), H1p(i), H1pGamma(i), maxpGamma(i), dof(i)] = FractureDarcyP0P0RT0_P0DG_2(eipsilon, i);
end
error0Rate = zeros(n,1); errorbRate = zeros(n,1); H1pRate = zeros(n,1);
H1pGammaRate = zeros(n,1); maxpGammaRate = zeros(n,1); digits(3);
for i = 2:n
    error0Rate(i) = vpa( log(error0(i-1)/error0(i))/log(2));
    errorbRate(i) = vpa( log(errorb(i-1)/errorb(i))/log(2));
    H1pRate(i) = vpa( log(H1p(i-1)/H1p(i))/log(2));
    H1pGammaRate(i) = vpa( log(H1pGamma(i-1)/H1pGamma(i))/log(2));
    maxpGammaRate(i) = vpa( log(maxpGamma(i-1)/maxpGamma(i))/log(2));
end
display('Table 3: error when eipsilon is 0.001');
colname = {'#DOF',         'error0',   'errorb', 'H1p',  'H1pGamma', 'maxpGamma'};
disptable(colname,dof(3:n),[], error0(3:n),'%0.4e', errorb(3:n),'%0.4e', H1p(3:n), '%0.4e',  H1pGamma(3:n),'%0.4e', maxpGamma(3:n),'%0.4e');
rate = [sum(error0Rate(3:n))/(n-2), sum(errorbRate(3:n))/(n-2), sum(H1pRate(3:n))/(n-2), sum(H1pGammaRate(3:n))/(n-2), sum(maxpGammaRate(3:n))/(n-2)];
for i = 1:5, if i<5, fprintf(' %0.2f  ', rate(i)); else fprintf(' %0.2f  \n', rate(i)); end; end
clear all