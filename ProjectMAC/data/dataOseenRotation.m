function [data, uh, vh, ph, f1h, f2h, gh, uI, vI, pI, width, nu] = dataOseenRotation(n, nu)

%% Exact data
data = struct('f1',@f1, 'f2',@f2, 'uexact',@uexact, 'vexact',@vexact, 'pexact',@pexact, 'omega',@omega);
kappa = 8; t = 3; s = 0.1;

%% Exact data
%     function r = uexact(x,y) 
%         r = s*exp(s*y)/(2*pi*(exp(s)-1)).*sin(2*pi*(exp(s*y)-1)/(exp(s)-1)).* ...
%             (1 - cos(2*pi*(exp(t*x)-1)/(exp(t)-1)));
%     end
%     function r = vexact(x,y)
%         r =  - t*exp(t*x)/(2*pi*(exp(t)-1)).*sin(2*pi*(exp(t*x)-1)/(exp(t)-1)).* ...
%              (1 - cos(2*pi*(exp(s*y)-1)/(exp(s)-1)));
%     end
%     function r = pexact(x,y)
%         r = t*s*exp(t*x).*exp(s*y)/((exp(t)-1)*(exp(s)-1)).* ...
%             sin(2*pi*(exp(t*x)-1)/(exp(t)-1)).* ...
%             sin(2*pi*(exp(s*y)-1)/(exp(s)-1));
%     end
%     function r = Pexact(x,y)
%         pp = pexact(x,y); uu = uexact(x,y); vv = vexact(x,y);
%         r = pp + 1/2 * (uu.*uu + vv.*vv);
%     end
%     function r = f1(x,y)
%         r = nu*((3*s^3*exp(2*s*y).*cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*(cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)) - 1))/(exp(s) - 1)^2 + ...
%             (s^3*exp(s*y).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*(cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)) - 1))/(2*pi*(exp(s) - 1)) - ...
%             (2*pi*s^3*exp(3*s*y).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*(cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)) - 1))/(exp(s) - 1)^3) - ...
%             nu*((s*t^2*exp(s*y).*exp(t*x).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)))/((exp(s) - 1)*(exp(t) - 1)) + ...
%             (2*pi*s*t^2*exp(s*y).*exp(2*t*x).*cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)))/((exp(s) - 1)*(exp(t) - 1)^2)) + ...
%             (t^3*exp(2*t*x).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).^2.*(cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)) - 1).^2)/(4*pi^2*(exp(t) - 1)^2) + ...
%             (t^3*exp(3*t*x).*cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*(cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)) - 1).^2)/(2*pi*(exp(t) - 1)^3) + ...
%             (s*t^2*exp(s*y + t*x).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)))/((exp(s) - 1)*(exp(t) - 1)) - ...
%             (t*exp(t*x).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*(cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)) - 1).*((s^2*exp(2*s*y).*cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*(cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)) - 1))/(exp(s) - 1)^2 + ...
%             (t^2*exp(2*t*x).*cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*(cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)) - 1))/(exp(t) - 1)^2 + ...
%             (s^2*exp(s*y).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*(cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)) - 1))/(2*pi*(exp(s) - 1)) + ...
%             (t^2*exp(t*x).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*(cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)) - 1))/(2*pi*(exp(t) - 1))))/(2*pi*(exp(t) - 1)) + ...
%             (2*pi*s*t^2*exp(s*y + t*x).*exp(t*x).*cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)))/((exp(s) - 1)*(exp(t) - 1)^2) - ...
%             (s^2*t*exp(t*x).*exp(2*s*y).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).^2.*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*(cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)) - 1))/(2*pi*(exp(s) - 1)^2*(exp(t) - 1));
%     end
%     function r = f2(x,y)
%         r = nu*((s^2*t*exp(s*y).*exp(t*x).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)))/((exp(s) - 1)*(exp(t) - 1)) + ...
%             (2*pi*s^2*t*exp(t*x).*exp(2*s*y).*cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)))/((exp(s) - 1)^2*(exp(t) - 1))) - ...
%             nu*((3*t^3*exp(2*t*x).*cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*(cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)) - 1))/(exp(t) - 1)^2 + ...
%             (t^3*exp(t*x).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*(cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)) - 1))/(2*pi*(exp(t) - 1)) - ...
%             (2*pi*t^3*exp(3*t*x).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*(cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)) - 1))/(exp(t) - 1)^3) + ...
%             (s^3*exp(2*s*y).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).^2.*(cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)) - 1).^2)/(4*pi^2*(exp(s) - 1)^2) + ...
%             (s^3*exp(3*s*y).*cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*(cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)) - 1).^2)/(2*pi*(exp(s) - 1)^3) + ...
%             (s^2*t*exp(s*y + t*x).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)))/((exp(s) - 1)*(exp(t) - 1)) - ...
%             (s*exp(s*y).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*(cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)) - 1).*((s^2*exp(2*s*y).*cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*(cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)) - 1))/(exp(s) - 1)^2 + ...
%             (t^2*exp(2*t*x).*cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*(cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)) - 1))/(exp(t) - 1)^2 + ...
%             (s^2*exp(s*y).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*(cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)) - 1))/(2*pi*(exp(s) - 1)) + ...
%             (t^2*exp(t*x).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*(cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)) - 1))/(2*pi*(exp(t) - 1))))/(2*pi*(exp(s) - 1)) + ...
%             (2*pi*s^2*t*exp(s*y + t*x).*exp(s*y).*cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)))/((exp(s) - 1)^2.*(exp(t) - 1)) - ...
%             (s*t^2*exp(s*y).*exp(2*t*x).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).^2.*(cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)) - 1))/(2*pi*(exp(s) - 1)*(exp(t) - 1)^2);
%     end
%     function r = omega(x,y)
%         r = (s^2*exp(2*s*y).*cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*(cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)) - 1))/(exp(s) - 1)^2 + ...
%             (t^2*exp(2*t*x).*cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*(cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)) - 1))/(exp(t) - 1)^2 + ...
%             (s^2*exp(s*y).*sin((2*pi*(exp(s*y) - 1))/(exp(s) - 1)).*(cos((2*pi*(exp(t*x) - 1))/(exp(t) - 1)) - 1))/(2*pi*(exp(s) - 1)) + ...
%             (t^2*exp(t*x).*sin((2*pi*(exp(t*x) - 1))/(exp(t) - 1)).*(cos((2*pi*(exp(s*y) - 1))/(exp(s) - 1)) - 1))/(2*pi*(exp(t) - 1));
%     end

    function r = uexact(x,y) 
        r = cos(2*pi*x).*sin(2*pi*y);
    end
    function r = vexact(x,y)
        r =  - sin(2*pi*x).*cos(2*pi*y);
    end
    function r = pexact(x,y)
        r = 0*(x+y);
    end
    function r = Pexact(x,y)
        pp = pexact(x,y); uu = uexact(x,y); vv = vexact(x,y);
        r = pp + 1/2 * (uu.*uu + vv.*vv);
%         r = pp;
    end
    function r = f1(x,y)
        u = uexact(x,y); v = vexact(x,y); dpx = 0; w = omega(x,y);
        Lapu = -8*pi^2*cos(2*pi*x).*sin(2*pi*y);
        dux = - 2*pi*sin(2*pi*x).*sin(2*pi*y);
        duy = 2*pi*cos(2*pi*x).*cos(2*pi*y);
        r = - nu*Lapu + u.*dux + v.*duy + dpx;
%         r = - nu*Lapu - w.*v + dpx;
    end
    function r = f2(x,y)
        u = uexact(x,y); v = vexact(x,y); dpy = 0; w = omega(x,y);
        Lapv = 8*pi^2*cos(2*pi*y).*sin(2*pi*x);
        dvx = - 2*pi*cos(2*pi*x).*cos(2*pi*y);
        dvy = 2*pi*sin(2*pi*x).*sin(2*pi*y);
        r = - nu*Lapv + u.*dvx + v.*dvy + dpy;
%         r = - nu*Lapv + w.*u + dpy;
    end
    function r = omega(x,y)
        r = - 4*pi*cos(2*pi*x).*cos(2*pi*y);
    end

%% Mesh
lef = 0; rig = 1; bot = 0; top = 1; h = (rig - lef) / n; width = rig - lef;
uh = zeros(n,n+1);  vh = zeros(n+1,n); ph = zeros(n,n);

[ux,uy] = meshgrid(lef: h: rig, top-h/2: -h: bot+h/2);
[vx,vy] = meshgrid(lef+h/2: h: rig-h/2, top: -h: bot);
[px,py] = meshgrid(lef+h/2: h: rig-h/2, top-h/2: -h: bot+h/2);

%% Impose Dirichlet for uh and vh
uLef = uexact(lef,uy(:,1)); uRig = uexact(rig,uy(:,end));
uh(:,[1 end]) = [uLef  uRig];
vTop = vexact(vx(1,:),top); vBot = vexact(vx(end,:),bot);
vh([1 end],:) = [vTop; vBot];

xx = lef: h: rig; j = 1:n+1; 
uTop(1,j) = uexact(xx(j),top); uBot(1,j) = uexact(xx(j),bot);
yy = top: -h: bot; i = 1:n+1; 
vLef(i,1) = vexact(lef,yy(i)); vRig(i,1) = vexact(rig,yy(i));

f1h = f1(ux(:),uy(:)); f1h = reshape(f1h,n,n+1);
f2h = f2(vx(:),vy(:)); f2h = reshape(f2h,n+1,n); gh = zeros(n,n);

%% Modify f1h and f2h for boundary condition for -Lap operator
f1h(:,[1 end]) = [uLef  uRig]; f2h([1 end],:) = [vTop;  vBot]; 

f1h(:,[2 end-1]) = f1h(:,[2 end-1]) + nu/h^2*[uLef  uRig]; 

f2h([2 end-1],:) = f2h([2 end-1],:) + nu/h^2*[vTop; vBot]; 

j = 2:n; f1h(1,j) = f1h(1,j) + nu/h^2*2*uTop(1,j); f1h(end,j) = f1h(end,j) + nu/h^2*2*uBot(1,j);

i = 2:n; f2h(i,1) = f2h(i,1) + nu/h^2*2*vLef(i,1); f2h(i,end) = f2h(i,end) + nu/h^2*2*vRig(i,1);


%% Modify f1h and f2h for boundary condition for omega*velocity operator
i = 1; j = 2:n; 
f1h(i,j) = f1h(i,j) + omega(ux(i,j),uy(i,j)).*vexact(vx(i,j-1),vy(i,j-1))/4;
f1h(i,j) = f1h(i,j) + omega(ux(i,j),uy(i,j)).*vexact(vx(i,j),vy(i,j))/4;

i = n; j = 2:n;
f1h(i,j) = f1h(i,j) + omega(ux(i,j),uy(i,j)).*vexact(vx(i+1,j-1),vy(i+1,j-1))/4;
f1h(i,j) = f1h(i,j) + omega(ux(i,j),uy(i,j)).*vexact(vx(i+1,j),vy(i+1,j))/4;

i = 2:n; j = 1; 
f2h(i,j) = f2h(i,j) - omega(vx(i,j),vy(i,j)).*uexact(ux(i-1,j),uy(i-1,j))/4;
f2h(i,j) = f2h(i,j) - omega(vx(i,j),vy(i,j)).*uexact(ux(i,j),uy(i,j))/4;

i = 2:n; j = n; 
f2h(i,j) = f2h(i,j) - omega(vx(i,j),vy(i,j)).*uexact(ux(i-1,j+1),uy(i-1,j+1))/4;
f2h(i,j) = f2h(i,j) - omega(vx(i,j),vy(i,j)).*uexact(ux(i,j+1),uy(i,j+1))/4;

%% Modify gh for boundary condition for -Div operator
gh(:,[1 end]) = gh(:,[1 end]) + [-uLef   uRig]/h; 
gh([1 end],:) = gh([1 end],:) + [vTop;  -vBot]/h; 

uI = uexact(ux(:),uy(:)); vI = vexact(vx(:),vy(:)); pI = Pexact(px(:),py(:)); 
  
clear ux uy vx vy px py xx yy

end