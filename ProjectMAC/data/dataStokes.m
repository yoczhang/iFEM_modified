function [uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, uI, vI, pI, h, width] = dataStokes(n)

%% Setup
lef = 0; rig = 1; top = 1; bot = 0; width = rig - lef;
h = (rig - lef) / n;
uh = zeros(n,n+1); vh = zeros(n+1,n); ph = zeros(n,n);

[ux,uy] = meshgrid(lef: h: rig, top-h/2: -h: bot+h/2);
[vx,vy] = meshgrid(lef+h/2: h: rig-h/2, top: -h: bot);
[px,py] = meshgrid(lef+h/2: h: rig-h/2, top-h/2: -h: bot+h/2);

yy = top-h/2: -h: bot+h/2; i = 1:n;
uh(i,1) = uexact(lef,yy(i)); uh(i,end) = uexact(rig,yy(i));

xx = lef+h/2: h: rig-h/2; j = 1:n;
vh(1,j) = vexact(xx(j),top); vh(end,j) = vexact(xx(j),bot);

% f1h = zeros(n, n+1); f2h = zeros(n+1, n); gh = zeros(n,n);
f1h = f1(ux(:),uy(:)); f1h = reshape(f1h,n,n+1);
f2h = f2(vx(:),vy(:)); f2h = reshape(f2h,n+1,n);
gh = zeros(n,n);

xx = lef: h: rig; i = 1:n+1;
uTop(1,i) = uexact(xx(i),top); uBot(1,i) = uexact(xx(i),bot);

yy = top: -h: bot; j = 1:n+1;
vLef(j,1) = vexact(lef,yy(j)); vRig(j,1) = vexact(rig,yy(j));

%% Interpolation
uI = uexact(ux(:),uy(:)); vI = vexact(vx(:),vy(:)); pI = pexact(px(:),py(:));


%% Exact solution
    function r = uexact(x,y) 
%         r = 20*x.*y.^3;
        r = x.^2.*(x-1).^2.*y.*(y-1).*(2.*y-1);
    end
    function r = vexact(x,y)
%         r = 5*x.^4-5*y.^4;
        r = -x.*(x-1).*(2.*x-1).*y.^2.*(y-1).^2;
    end
    function r = pexact(x,y)
%         r = 60*x.^2.*y - 20*y.^3;
        r = (2.*x-1).*(2.*y-1);
    end
    function r = f1(x,y)
%         r = zeros(size(x,1),1);
        nu = 1; r = 2*(2*y - 1).*(- 3*nu*x.^4 + 6*nu*x.^3 - 6*nu*x.^2.*y.^2 + 6*nu*x.^2.*y - 3*nu*x.^2 + 6*nu*x.*y.^2 - 6*nu*x.*y - nu*y.^2 + nu*y + 1);
    end
    function r = f2(x,y)
%         r = zeros(size(x,1),1);
        nu = 1; r = 2*(2*x - 1).*(6*nu*x.^2.*y.^2 - 6*nu*x.^2.*y + nu*x.^2 - 6*nu*x.*y.^2 + 6*nu*x.*y - nu*x + 3*nu*y.^4 - 6*nu*y.^3 + 3*nu*y.^2 + 1);
    end

end