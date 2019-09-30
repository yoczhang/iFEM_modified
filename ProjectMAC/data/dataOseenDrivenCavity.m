function  [data, uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig] = dataOseenDrivenCavity(n)

%% Exact data
data = struct('flowa',@flowa, 'flowb',@flowb);

    function r = flowa(x,y)
        r = 8*x.*(x-1).*(1-2*y);
    end
    function r = flowb(x,y)
        r = 8*(2*x-1).*y.*(y-1);
    end


%% Discrete data
uh = zeros(n,n+1); vh = zeros(n+1,n); ph = zeros(n,n);
f1h = zeros(n, n+1); f2h = zeros(n+1, n); gh = zeros(n,n);
uTop = ones(1,n+1); uBot = zeros(1,n+1); vLef = zeros(n+1,1); vRig = zeros(n+1,1);

end