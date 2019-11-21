function [uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, uI, vI, pI, h, width] = gmres_dataSaddlePoint(n, mu, gamma)
%
%

    %% u
    function r = uexact(x,y) 
        r = (1 - cos(2*pi*x)).*sin(2*pi*y);
    end
    function r = uxexact(x,y) 
        r = 2*pi*sin(2*pi*x).*sin(2*pi*y);
    end
    function r = uxxexact(x,y) 
        r = 4*pi^2*cos(2*pi*x).*sin(2*pi*y);
    end
    function r = uyexact(x,y) 
        r = -2*pi*cos(2*pi*y).*(cos(2*pi*x) - 1);
    end
    function r = uyyexact(x,y) 
        r = 4*pi^2*sin(2*pi*y).*(cos(2*pi*x) - 1);
    end
    
    %% v
    function r = vexact(x,y)
        r = (cos(2*pi*y) - 1).*sin(2*pi*x);
    end
    function r = vxexact(x,y)
        r = 2*pi*cos(2*pi*x).*(cos(2*pi*y) - 1);
    end
    function r = vxxexact(x,y)
        r = -4*pi^2*sin(2*pi*x).*(cos(2*pi*y) - 1);
    end
    function r = vyexact(x,y)
        r = -2*pi*sin(2*pi*x).*sin(2*pi*y);
    end
    function r = vyyexact(x,y)
        r = -4*pi^2*cos(2*pi*y).*sin(2*pi*x);
    end

    %% p
    function r = pexact(x,y)
        r = x.^3 / 3 - 1 / 12;
    end
    function r = pxexact(x,y)
        r = x.^2;
    end
    function r = pyexact(x,y)
        r = 0.*x;
    end

    %% f
    function r = f1(x,y)
        u_ = uexact(x,y);
        uxx_ = uxxexact(x,y);
        uyy_ = uyyexact(x,y);
        px_ = pxexact(x,y);
        r = -mu*(uxx_+uyy_) + gamma*u_ + px_;
    end
    function r = f2(x,y)
        v_ = vexact(x,y);
        vxx_ = vxxexact(x,y);
        vyy_ = vyyexact(x,y);
        py_ = pyexact(x,y);
        r = -mu*(vxx_+vyy_) + gamma*v_ + py_;
    end
    function r = g(x,y)
        ux_ = uxexact(x,y);
        vy_ = vyexact(x,y);
        p_ = pexact(x,y);
        r = -(ux_ + vy_ ) - p_;
    end


%% Discrete data
lef = 0; rig = 1; top = 1; bot = 0; h = (rig - lef) / n; width = rig - lef;
uh = zeros(n,n+1);   vh = zeros(n+1,n); ph = zeros(n,n);

[ux,uy] = meshgrid(lef: h: rig, top-h/2: -h: bot+h/2);
[vx,vy] = meshgrid(lef+h/2: h: rig-h/2, top: -h: bot);
[px,py] = meshgrid(lef+h/2: h: rig-h/2, top-h/2: -h: bot+h/2);

yy = top-h/2: -h: bot+h/2; i = 1:n; uh(i,1) = uexact(lef,yy(i)); uh(i,end) = uexact(rig,yy(i));

xx = lef+h/2: h: rig-h/2; j = 1:n; vh(1,j) = vexact(xx(j),top); vh(end,j) = vexact(xx(j),bot);

xx = lef: h: rig; i = 1:n+1; uTop(1,i) = uexact(xx(i),top); uBot(1,i) = uexact(xx(i),bot);

yy = top: -h: bot; j = 1:n+1; vLef(j,1) = vexact(lef,yy(j)); vRig(j,1) = vexact(rig,yy(j));

uI = uexact(ux(:),uy(:)); vI = vexact(vx(:),vy(:)); pI = pexact(px(:),py(:)); 

%- generate the random 0,1 matrix
% try 
%     saveFilename = ['kappa_k_',num2str(n)];
%     load(saveFilename, 'kappa_k_u', 'kappa_k_v')
% catch
%     randmat1 = rand(n*(n+1), 1);
%     randmat2 = zeros(n*(n+1),1);
%     randmat2(randmat1<=0.3) = 1;
%     randmat2(randmat1>0.3) = 0;
% 
%     kappa_k_u = reshape(randmat2, n, n+1);
%     kappa_k_v = reshape(randmat2, n+1, n);
% end
randmat1 = rand(n*(n+1), 1);
randmat2 = zeros(n*(n+1),1);
randmat2(randmat1<=0.3) = 0;
randmat2(randmat1>0.3) = 0;

kappa_k_u = reshape(randmat2, n, n+1);
kappa_k_v = reshape(randmat2, n+1, n);
saveFilename = ['kappa_k_',num2str(n)];
save( saveFilename, 'kappa_k_u', 'kappa_k_v', 'mu', 'gamma', 'uTop', 'uBot', 'vLef', 'vRig')

f1h = f1(ux(:),uy(:)); 
f1h = reshape(f1h,n,n+1); 
uI = reshape(uI,n,n+1); 
f1h = f1h + kappa_k_u.*uI;

f2h = f2(vx(:),vy(:)); 
f2h = reshape(f2h,n+1,n);  
vI = reshape(vI,n+1,n); 
f2h = f2h + kappa_k_v.*vI;
gh = g(px(:),py(:)); gh = reshape(gh,n,n);


clear ux uy vx vy px py xx yy
end