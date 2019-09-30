function  [data, uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, uI, vI, pI] = dataOseen(n, mu)

%% Exact data
data = struct('f1',@f1, 'f2',@f2, 'uexact',@uexact, 'vexact',@vexact, 'pexact',@pexact, 'flowa',@flowa, 'flowb',@flowb);

    function r = uexact(x,y) 
        r = (1 - cos(2*pi*x)).*sin(2*pi*y);
    end
    function r = vexact(x,y)
        r = (cos(2*pi*y) - 1).*sin(2*pi*x);
    end
    function r = pexact(x,y)
        r = x.^3 / 3 - 1 / 12;
    end
    function r = f1(x,y)
        r = x.^2 - 4*pi^2*mu*sin(2*pi*y).*(cos(2*pi*x) - 1) - ...
            4*pi^2*mu*cos(2*pi*x).*sin(2*pi*y) + ...
            2*pi*x.*sin(2*pi*x).*sin(2*pi*y).^2 - ...
            2*pi*y.*cos(2*pi*y).*sin(2*pi*x).*(cos(2*pi*x) - 1);
    end
    function r = f2(x,y)
        r = 4*pi^2*mu*sin(2*pi*x).*(cos(2*pi*y) - 1) + ...
            4*pi^2*mu*cos(2*pi*y).*sin(2*pi*x) - ...
            2*pi*y.*sin(2*pi*x).^2.*sin(2*pi*y) + ...
            2*pi*x.*cos(2*pi*x).*sin(2*pi*y).*(cos(2*pi*y) - 1);
    end
    function r = flowa(x,y)
        r = x.*sin(2*pi*y);
    end
    function r = flowb(x,y)
        r = y.*sin(2*pi*x);
    end

%% Discrete data
lef = 0; rig = 1; top = 1; bot = 0; h = (rig - lef) / n; 
uh = zeros(n,n+1);   vh = zeros(n+1,n); ph = zeros(n,n);

[ux,uy] = meshgrid(lef: h: rig, top-h/2: -h: bot+h/2);
[vx,vy] = meshgrid(lef+h/2: h: rig-h/2, top: -h: bot);
[px,py] = meshgrid(lef+h/2: h: rig-h/2, top-h/2: -h: bot+h/2);

yy = top-h/2: -h: bot+h/2; i = 1:n; uh(i,1) = uexact(lef,yy(i)); uh(i,end) = uexact(rig,yy(i));

xx = lef+h/2: h: rig-h/2; j = 1:n; vh(1,j) = vexact(xx(j),top); vh(end,j) = vexact(xx(j),bot);

f1h = f1(ux(:),uy(:)); f1h = reshape(f1h,n,n+1);
f2h = f2(vx(:),vy(:)); f2h = reshape(f2h,n+1,n); gh = zeros(n,n);

xx = lef: h: rig; i = 1:n+1; uTop(1,i) = uexact(xx(i),top); uBot(1,i) = uexact(xx(i),bot);

yy = top: -h: bot; j = 1:n+1; vLef(j,1) = vexact(lef,yy(j)); vRig(j,1) = vexact(rig,yy(j));

uI = uexact(ux(:),uy(:)); vI = vexact(vx(:),vy(:)); pI = pexact(px(:),py(:)); 
        
clear ux uy vx vy px py xx yy
end