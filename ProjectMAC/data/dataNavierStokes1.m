function  [data, uh, vh, ph, f1h, f2h, gh, uLef, uRig, uTop, uBot, vLef, vRig, vTop, vBot, uI, vI, pI] = dataNavierStokes1(n, nu, initialType)

%% Exact data
data = struct('f1',@f1, 'f2',@f2, 'uexact',@uexact, 'vexact',@vexact, 'pexact',@pexact, 'emptyFunc',@emptyFunc, 'flowa',@flowa, 'flowb',@flowb);

    function r = uexact(x,y) 
%         r = (1 - cos(2*pi*x)).*sin(2*pi*y);
        r = 0*(x+y);
    end
    function r = vexact(x,y)
%         r = (cos(2*pi*y) - 1).*sin(2*pi*x);
        r = 0*(x+y);
    end
    function r = pexact(x,y)
%         r = x.^3 / 3 - 1 / 12;
        r = 0*(x+y);
    end
    function r = f1(x,y) 
%         r = x.^2 + 8*pi^2*nu*sin(pi*x).^2.*sin(2*pi*y) + ...
%             4*pi*sin(pi*x).^2.*sin(2*pi*x).*sin(2*pi*y).^2 + ...
%             4*pi^2*nu*sin(2*pi*y).*(2*sin(pi*x).^2 - 1) + ...
%             8*pi*sin(pi*x).^2.*sin(2*pi*x).*sin(pi*y).^2.*(2*sin(pi*y).^2 - 1);
        r = 0*(x+y);
    end
    function r = f2(x,y)
%         r = 4*pi*sin(2*pi*x).^2.*sin(pi*y).^2.*sin(2*pi*y) - ...
%             8*pi^2*nu*sin(2*pi*x).*sin(pi*y).^2 - ...
%             4*pi^2*nu*sin(2*pi*x).*(2*sin(pi*y).^2 - 1) + ...
%             8*pi*sin(pi*x).^2.*sin(pi*y).^2.*sin(2*pi*y).*(2*sin(pi*x).^2 - 1);
        r = 0*(x+y);
    end


%% Discrete data
lef = 0; rig = 1; top = 1; bot = 0; h = (rig - lef) / n; 
if strcmp(initialType,'zeroInitial')
    uh = zeros(n,n+1);   vh = zeros(n+1,n); ph = zeros(n,n);
elseif strcmp(initialType,'randomInitial')
    uh = rand(n,n+1); vh = rand(n+1,n); ph = rand(n,n); ph = ph - mean(ph(:));
end

[ux,uy] = meshgrid(lef: h: rig, top-h/2: -h: bot+h/2);
[vx,vy] = meshgrid(lef+h/2: h: rig-h/2, top: -h: bot);
[px,py] = meshgrid(lef+h/2: h: rig-h/2, top-h/2: -h: bot+h/2);

%% Impose Dirichlet for uh and vh
uLef = uexact(lef,uy(:,1)); uRig = uexact(rig,uy(:,end)); uh(:,[1 end]) = [uLef  uRig];
vTop = vexact(vx(1,:),top); vBot = vexact(vx(end,:),bot); vh([1 end],:) = [vTop; vBot];

xx = lef: h: rig; j = 1:n+1;  uTop(1,j) = uexact(xx(j),top); uBot(1,j) = uexact(xx(j),bot);
yy = top: -h: bot; i = 1:n+1; vLef(i,1) = vexact(lef,yy(i)); vRig(i,1) = vexact(rig,yy(i));

f1h = f1(ux(:),uy(:)); f1h = reshape(f1h,n,n+1);
f2h = f2(vx(:),vy(:)); f2h = reshape(f2h,n+1,n); gh = zeros(n,n);

%% Modify f1h and f2h for boundary condition for -Lap operator
% f1h(:,[1 end]) = [uLef  uRig]; f2h([1 end],:) = [vTop;  vBot]; 
% 
% f1h(:,[2 end-1]) = f1h(:,[2 end-1]) + nu/h^2*[uLef  uRig]; 
% f2h([2 end-1],:) = f2h([2 end-1],:) + nu/h^2*[vTop; vBot]; 
% 
% j = 2:n; f1h(1,j) = f1h(1,j) + nu/h^2*2*uTop(1,j); f1h(end,j) = f1h(end,j) + nu/h^2*2*uBot(1,j);
% i = 2:n; f2h(i,1) = f2h(i,1) + nu/h^2*2*vLef(i,1); f2h(i,end) = f2h(i,end) + nu/h^2*2*vRig(i,1);

%% Modify gh for boundary condition for -Div operator
gh(:,[1 end]) = gh(:,[1 end]) + [-uLef   uRig]/h; 
gh([1 end],:) = gh([1 end],:) + [vTop;  -vBot]/h; 
gh = gh - mean(gh(:));

uI = uexact(ux(:),uy(:)); vI = vexact(vx(:),vy(:)); pI = pexact(px(:),py(:)); 
  
clear ux uy vx vy px py xx yy
end