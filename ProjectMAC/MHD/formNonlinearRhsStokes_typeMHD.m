function [f1, f2, g1, g2] = formNonlinearRhsStokes_typeMHD(data, n, h, u1, u2, B1, B2, t, lef, rig, bot, top)

f1 = zeros(n,n+1); f2 = zeros(n+1,n); g1 = zeros(n,n+1); g2 = zeros(n+1,n);
x = lef: h: rig;   y = (top: -h: bot)'; 

%% Compute J = dB2/dx - dB1/dy
J = zeros(n+1,n+1);
i = 2:n; j = 2:n; J(i,j) = (B2(i,j) - B2(i,j-1))/h - (B1(i-1,j) - B1(i,j))/h;

J(i,1) = (2*B2(i,1) - 2*data.B2exact(lef,y(2:n),t))/h - (B1(i-1,1) - B1(i,1))/h;
J(i,n+1) = (2*data.B2exact(rig,y(2:n),t) - B2(i,n))/h - (B1(i-1,n+1) - B1(i,n+1))/h;

J(1,j) = (B2(1,j) - B2(1,j-1))/h - (2*data.B1exact(x(2:n),top,t) - 2*B1(1,j))/h;
J(n+1,j) = (B2(n+1,j) - B2(n+1,j-1))/h - (2*B1(n,j) - 2*data.B1exact(x(2:n),bot,t))/h;

J(1,1) = (2*B2(1,1) - data.B2exact(lef,top,t))/h - (2*data.B1exact(lef,top,t) - 2*B1(1,1))/h;
J(1,n+1) = (2*data.B2exact(rig,top,t) - 2*B2(1,n))/h - (2*data.B1exact(rig,top,t) - 2*B1(1,n+1))/h;
J(n+1,1) = (2*B2(n+1,1) - 2*data.B2exact(lef,bot,t))/h - (2*B1(n,1) - 2*data.B1exact(lef,bot,t))/h;
J(n+1,n+1) = (2*data.B2exact(rig,bot,t) - 2*B2(n+1,n))/h - (2*B1(n,n+1) - 2*data.B1exact(rig,bot,t))/h;

%% Compute du1/dx, du1/dy, du2/dy, dB1/dy, dB2/dy at centers of vertical edges
du1dxVer = zeros(n,n+1); du1dyVer = zeros(n,n+1); du2dyVer = zeros(n,n+1); dB1dyVer = zeros(n,n+1); dB2dyVer = zeros(n,n+1);

i = 1:n; j = 2:n; du1dxVer(i,j) = (u1(i,j+1) - u1(i,j-1)) / (2*h);

i = 2:n-1; j = 2:n; du1dyVer(i,j) = (u1(i-1,j) - u1(i+1,j)) / (2*h);
du1dyVer(1,j) = (2*data.u1exact(x(2:n),top,t) - u1(1,j) - u1(2,j)) / (2*h);
du1dyVer(n,j) = (u1(n-1,j) + u1(n,j) - 2*data.u1exact(x(2:n),bot,t)) / (2*h);

i = 2:n-1; j = 2:n; dB1dyVer(i,j) = (B1(i-1,j) - B1(i+1,j)) / (2*h);
dB1dyVer(1,j) = (2*data.B1exact(x(2:n),top,t) - B1(1,j) - B1(2,j)) / (2*h);
dB1dyVer(n,j) = (B1(n-1,j) + B1(n,j) - 2*data.B1exact(x(2:n),bot,t)) / (2*h);

i = 2:n-1; j = 2:n;
du2dyVer(i,j) = ((u2(i-1,j-1)+u2(i,j-1)+u2(i-1,j)+u2(i,j))/4 - (u2(i+1,j-1)+u2(i+2,j-1)+u2(i+1,j)+u2(i+2,j))/4) / (2*h);
du2dyVer(1,j) = (2*data.u2exact(x(2:n),top,t) - (u2(1,j-1)+u2(2,j-1)+u2(1,j)+u2(2,j))/4 - (u2(2,j-1)+u2(3,j-1)+u2(2,j)+u2(3,j))/4) / (2*h);
du2dyVer(n,j) = ((u2(n-1,j-1)+u2(n,j-1)+u2(n-1,j)+u2(n,j))/4 + (u2(n,j-1)+u2(n+1,j-1)+u2(n,j)+u2(n+1,j))/4 - 2*data.u2exact(x(2:n),bot,t)) / (2*h);

i = 2:n-1; j = 2:n;
dB2dyVer(i,j) = ((B2(i-1,j-1)+B2(i,j-1)+B2(i-1,j)+B2(i,j))/4 - (B2(i+1,j-1)+B2(i+2,j-1)+B2(i+1,j)+B2(i+2,j))/4) / (2*h);
dB2dyVer(1,j) = (2*data.B2exact(x(2:n),top,t) - (B2(1,j-1)+B2(2,j-1)+B2(1,j)+B2(2,j))/4 - (B2(2,j-1)+B2(3,j-1)+B2(2,j)+B2(3,j))/4) / (2*h);
dB2dyVer(n,j) = ((B2(n-1,j-1)+B2(n,j-1)+B2(n-1,j)+B2(n,j))/4 + (B2(n,j-1)+B2(n+1,j-1)+B2(n,j)+B2(n+1,j))/4 - 2*data.B2exact(x(2:n),bot,t)) / (2*h);

%% Compute du1/dx, du2/dx, du2/dy, dB1/dx, dB2/dx at centers of horizontal edges
du1dxHor = zeros(n+1,n); du2dxHor = zeros(n+1,n); du2dyHor = zeros(n+1,n); dB1dxHor = zeros(n+1,n); dB2dxHor = zeros(n+1,n);

i = 2:n; j = 2:n-1;
du1dxHor(i,j) = ((u1(i-1,j-1)+u1(i,j-1)+u1(i-1,j)+u1(i,j))/4 - (u1(i-1,j+1)+u1(i,j+1)+u1(i-1,j+2)+u1(i,j+2))/4) / (2*h);
du1dxHor(i,1) = ((u1(i-1,2)+u1(i,2)+u1(i-1,3)+u1(i,3))/4 + (u1(i-1,1)+u1(i,1)+u1(i-1,2)+u1(i,2))/4 - 2*data.u1exact(lef,y(2:n),t)) / (2*h);
du1dxHor(i,n) = (2*data.u1exact(rig,y(2:n),t) - (u1(i-1,n-1)+u1(i,n-1)+u1(i-1,n)+u1(i,n))/4 - (u1(i-1,n)+u1(i,n)+u1(i-1,n+1)+u1(i,n+1))/4) / (2*h);

i = 2:n; j = 2:n-1;
dB1dxHor(i,j) = ((B1(i-1,j-1)+B1(i,j-1)+B1(i-1,j)+B1(i,j))/4 - (B1(i-1,j+1)+B1(i,j+1)+B1(i-1,j+2)+B1(i,j+2))/4) / (2*h);
dB1dxHor(i,1) = ((B1(i-1,2)+B1(i,2)+B1(i-1,3)+B1(i,3))/4 + (B1(i-1,1)+B1(i,1)+B1(i-1,2)+B1(i,2))/4 - 2*data.B1exact(lef,y(2:n),t)) / (2*h);
dB1dxHor(i,n) = (2*data.B1exact(rig,y(2:n),t) - (B1(i-1,n-1)+B1(i,n-1)+B1(i-1,n)+B1(i,n))/4 - (B1(i-1,n)+B1(i,n)+B1(i-1,n+1)+B1(i,n+1))/4) / (2*h);

i = 2:n; j = 2:n-1;
du2dxHor(i,j) = (u2(i,j+1) - u2(i,j-1)) / (2*h);
du2dxHor(i,1) = (u2(i,2) - (2*data.u2exact(lef,y(2:n),t) - u2(i,1))) / (2*h);
du2dxHor(i,n) = ((2*data.u2exact(rig,y(2:n),t) - u2(i,n)) - u2(i,n-1)) / (2*h);

i = 2:n; j = 2:n-1;
dB2dxHor(i,j) = (B2(i,j+1) - B2(i,j-1)) / (2*h);
dB2dxHor(i,1) = (B2(i,2) - (2*data.B2exact(lef,y(2:n),t) - B2(i,1))) / (2*h);
dB2dxHor(i,n) = ((2*data.B2exact(rig,y(2:n),t) - B2(i,n)) - B2(i,n-1)) / (2*h);

i = 2:n; j = 1:n;
du2dyHor(i,j) = (u2(i-1,j) - u2(i+1,j)) / (2*h);

%% Lastly, to compute f1, f2, g1, g2
i = 1:n; j = 2:n;
f1(i,j) = - u1(i,j) .* du1dxVer(i,j) - (u2(i,j-1)+u2(i+1,j-1)+u2(i,j)+u2(i+1,j))/4 .* du1dyVer(i,j) - ...
          (J(i,j)+J(i+1,j))/2 .* (B2(i,j-1)+B2(i+1,j-1)+B2(i,j)+B2(i+1,j))/4;
g1(i,j) = du1dyVer(i,j) .* (B2(i,j-1)+B2(i+1,j-1)+B2(i,j)+B2(i+1,j))/4 + ...
          u1(i,j) .* dB2dyVer(i,j) - ...
          du2dyVer(i,j) .* B1(i,j) - ...
          (u2(i,j-1)+u2(i+1,j-1)+u2(i,j)+u2(i+1,j))/4 .* dB1dyVer(i,j);
      
i = 2:n; j = 1:n;
f2(i,j) = - (u1(i-1,j)+u1(i,j)+u1(i-1,j+1)+u1(i,j+1))/4 .* du2dxHor(i,j) - u2(i,j) .* du2dyHor(i,j) + ...
          (J(i,j)+J(i,j+1))/2 .* (B1(i-1,j)+B1(i,j)+B1(i-1,j+1)+B1(i,j+1))/4;
g2(i,j) = - du1dxHor(i,j) .* B2(i,j) - ...
          (u1(i-1,j)+u1(i,j)+u1(i-1,j+1)+u1(i,j+1))/4 .* dB2dxHor(i,j) + ...
          du2dxHor(i,j) .* (B1(i-1,j)+B1(i,j)+B1(i-1,j+1)+B1(i,j+1))/4 + ...
          u2(i,j) .* dB1dxHor(i,j);
end