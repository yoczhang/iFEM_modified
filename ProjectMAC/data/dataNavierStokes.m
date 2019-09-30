function  [data, uh, vh, ph, f1h, f2h, gh, uLef, uRig, uTop, uBot, vLef, vRig, vTop, vBot, uI, vI, pI] = dataNavierStokes(n, viscosity)

%% Exact data
data = struct('f1',@f1, 'f2',@f2, 'uexact',@uexact, 'vexact',@vexact, 'pexact',@pexact, 'emptyFunc',@emptyFunc, 'flowa',@flowa, 'flowb',@flowb);

    function r = uexact(x,y) 
        r = 10*x.^2.*(x-1).^2.*y.*(y-1).*(2*y-1);
%         r = (1 - cos(2*pi*x)).*sin(2*pi*y);
%         r = (1 - cos(2*pi*x)).*sin(2*pi*y);
%         r = (y==1)*1 + (y<0)*0;
%          r = cos(2*pi*x).*sin(2*pi*y);
    end
    function r = vexact(x,y)
        r = -10*x.*(x-1).*(2*x-1).*y.^2.*(y-1).^2;
%         r = (cos(2*pi*y) - 1).*sin(2*pi*x);
%         r = (cos(2*pi*y) - 1).*sin(2*pi*x);
%         r = 0*(x+y);
%          r = - sin(2*pi*x).*cos(2*pi*y);
    end
    function r = pexact(x,y)
        r = 10*(2*x-1).*(2*y-1);
%         r = x.^3 / 3 - 1 / 12;
%         r = x.^3 / 3 - 1 / 12;
%         r = 0*(x+y);
%          r = 0*(x+y);
    end
    function r = f1(x,y)
        r = 40*y - 60*viscosity*x.^2.*(2*y - 1).*(x - 1).^2 - 20*viscosity*y.*(6*x.^2 - 6*x + 1).*(2*y.^2 - 3*y + 1) - ...
            100*x.^3.*y.^2.*(2*x - 1).*(x - 1).^3.*(y - 1).^2.*(6*y.^2 - 6*y + 1) + ...
            200*x.^3.*y.^2.*(2*y - 1).*(x - 1).^2.*(y - 1).*(2*x.^2 - 3*x + 1).*(2*y.^2 - 3*y + 1) - 20;
%         r = x.^2 + 8*pi^2*viscosity*sin(pi*x).^2.*sin(2*pi*y) + ...
%             4*pi*sin(pi*x).^2.*sin(2*pi*x).*sin(2*pi*y).^2 + ...
%             4*pi^2*viscosity*sin(2*pi*y).*(2*sin(pi*x).^2 - 1) + ...
%             8*pi*sin(pi*x).^2.*sin(2*pi*x).*sin(pi*y).^2.*(2*sin(pi*y).^2 - 1);
%         r = x.^2 - 4*pi^2*viscosity*sin(2*pi*y).*(cos(2*pi*x) - 1) - ...
%             4*pi^2*viscosity*cos(2*pi*x).*sin(2*pi*y) + ...
%             2*pi*x.*sin(2*pi*x).*sin(2*pi*y).^2 - ...
%             2*pi*y.*cos(2*pi*y).*sin(2*pi*x).*(cos(2*pi*x) - 1);
%         r = 0*(x+y);
%          r = -2*pi*cos(2*pi*x).*(sin(2*pi*x) - 4*pi*viscosity*sin(2*pi*y));
    end
    function r = f2(x,y)
        r = 40*x + 60*viscosity*y.^2.*(2*x - 1).*(y - 1).^2 + 20*viscosity*x.*(2*x.^2 - 3*x + 1).*(6*y.^2 - 6*y + 1) - ...
            100*x.^2.*y.^3.*(2*y - 1).*(x - 1).^2.*(y - 1).^3.*(6*x.^2 - 6*x + 1) + ...
            200*x.^2.*y.^3.*(2*x - 1).*(x - 1).*(y - 1).^2.*(2*x.^2 - 3*x + 1).*(2*y.^2 - 3*y + 1) - 20;
%         r = 4*pi*sin(2*pi*x).^2.*sin(pi*y).^2.*sin(2*pi*y) - ...
%             8*pi^2*viscosity*sin(2*pi*x).*sin(pi*y).^2 - ...
%             4*pi^2*viscosity*sin(2*pi*x).*(2*sin(pi*y).^2 - 1) + ...
%             8*pi*sin(pi*x).^2.*sin(pi*y).^2.*sin(2*pi*y).*(2*sin(pi*x).^2 - 1);
%         r = 4*pi^2*viscosity*sin(2*pi*x).*(cos(2*pi*y) - 1) + ...
%             4*pi^2*viscosity*cos(2*pi*y).*sin(2*pi*x) - ...
%             2*pi*y.*sin(2*pi*x).^2.*sin(2*pi*y) + ...
%             2*pi*x.*cos(2*pi*x).*sin(2*pi*y).*(cos(2*pi*y) - 1);
%         r = 0*(x+y);
%          r = -2*pi*cos(2*pi*y).*(sin(2*pi*y) + 4*pi*viscosity*sin(2*pi*x));
    end
    function r = emptyFunc(x,y)
        r = 0*(x+y);
    end
    function r = flowa(x,y)
        r = 10*x.^2.*(x-1).^2.*y.*(y-1).*(2*y-1);
%         r = (1 - cos(2*pi*x)).*sin(2*pi*y);
%         r = x.*sin(2*pi*y);
%         r = 8*x.*(x-1).*(1-2*y);
%          r = cos(2*pi*x).*sin(2*pi*y);
    end
    function r = flowb(x,y)
        r = -10*x.*(x-1).*(2*x-1).*y.^2.*(y-1).^2;
%         r = (cos(2*pi*y) - 1).*sin(2*pi*x);
%         r = y.*sin(2*pi*x);
%         r = 8*(2*x-1).*y.*(y-1);
%          r = - sin(2*pi*x).*cos(2*pi*y);
    end

%% Discrete data
lef = 0; rig = 1; top = 1; bot = 0; h = (rig - lef) / n; 
uh = zeros(n,n+1);   vh = zeros(n+1,n); ph = zeros(n,n);

[ux,uy] = meshgrid(lef: h: rig, top-h/2: -h: bot+h/2);
[vx,vy] = meshgrid(lef+h/2: h: rig-h/2, top: -h: bot);
[px,py] = meshgrid(lef+h/2: h: rig-h/2, top-h/2: -h: bot+h/2);

%% Impose Dirichlet for uh and vh
uLef = uexact(lef,uy(:,1)); uRig = uexact(rig,uy(:,end)); uh(:,[1 end]) = [uLef  uRig];
vTop = vexact(vx(1,:),top); vBot = vexact(vx(end,:),bot); vh([1 end],:) = [vTop; vBot];

xx = lef:  h: rig; j = 1:n+1; uTop(1,j) = uexact(xx(j),top); uBot(1,j) = uexact(xx(j),bot);
yy = top: -h: bot; i = 1:n+1; vLef(i,1) = vexact(lef,yy(i)); vRig(i,1) = vexact(rig,yy(i));

f1h = f1(ux(:),uy(:)); f1h = reshape(f1h,n,n+1);
f2h = f2(vx(:),vy(:)); f2h = reshape(f2h,n+1,n); gh = zeros(n,n);

%% Modify f1h and f2h for boundary condition for -Lap operator
% f1h(:,[1 end]) = [uLef  uRig]; f2h([1 end],:) = [vTop;  vBot]; 
% 
% f1h(:,[2 end-1]) = f1h(:,[2 end-1]) + viscosity/h^2*[uLef  uRig]; 
% f2h([2 end-1],:) = f2h([2 end-1],:) + viscosity/h^2*[vTop; vBot]; 
% 
% j = 2:n; f1h(1,j) = f1h(1,j) + viscosity/h^2*2*uTop(1,j); f1h(end,j) = f1h(end,j) + viscosity/h^2*2*uBot(1,j);
% i = 2:n; f2h(i,1) = f2h(i,1) + viscosity/h^2*2*vLef(i,1); f2h(i,end) = f2h(i,end) + viscosity/h^2*2*vRig(i,1);

%% Modify gh for boundary condition for -Div operator
gh(:,[1 end]) = gh(:,[1 end]) + [-uLef   uRig]/h; 
gh([1 end],:) = gh([1 end],:) + [vTop;  -vBot]/h; 

uI = uexact(ux(:),uy(:)); vI = vexact(vx(:),vy(:)); pI = pexact(px(:),py(:)); 
  
clear ux uy vx vy px py xx yy
end