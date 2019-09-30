function data = dataRotation2(nu)

%% Exact data
data = struct('f',@f, 'exactu',@exactu, 'exactp',@exactp, 'littlep',@littlep, 'omega',@omega, 'g_D',@g_D, 'diffCo',@diffCo);
kappa = 8;

%% Exact data
    function r = f(p)
        r = zeros(size(p,1),2);
    end 
    function r = omega(p)
        x = p(:,1); y = p(:,2);
        r = 2*kappa*x.*(x-1) + 2*kappa*y.*(y-1);
    end
    function r = g_D(p)
        x = p(:,1); y = p(:,2);
        r = zeros(size(p,1),2);
        bdtop = (abs(y-1)==0);
        r(bdtop,1) = 1;
    end
    function r = exactu(p)
        r = zeros(size(p,1),2);
    end
    function r = exactp(p)
        r = zeros(size(p,1),1);
    end
    function r = diffCo(p)
        r = nu*ones(size(p,1),1);
    end
end
