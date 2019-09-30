function formAnalyticalSolutions

syms x y nu t s

% u = 20*x.*y.^3;
% v = 5*x.^4-5*y.^4;
% P = 60*x.^2.*y - 20*y.^3;
% p = P - 1/2*(u.*u + v.*v);

u = s*exp(s*y)/(2*pi*(exp(s)-1)).*sin(2*pi*(exp(s*y)-1)/(exp(s)-1)).* ...
    (1 - cos(2*pi*(exp(t*x)-1)/(exp(t)-1)));
v = - t*exp(t*x)/(2*pi*(exp(t)-1)).*sin(2*pi*(exp(t*x)-1)/(exp(t)-1)).* ...
    (1 - cos(2*pi*(exp(s*y)-1)/(exp(s)-1)));
p = t*s*exp(t*x).*exp(s*y)/((exp(t)-1)*(exp(s)-1))* ...
    sin(2*pi*(exp(t*x)-1)/(exp(t)-1)).* ...
    sin(2*pi*(exp(s*y)-1)/(exp(s)-1));
P = simple(p + 1/2*(u.*u + v.*v));

w  = simple(diff(v,x) - diff(u,y));
P, w
divValue = simple(diff(u,x) + diff(v,y));
divValue
%% Rotation form
f1 = simple( - nu*diff(u,x,2) - nu*diff(u,y,2) - w*v + diff(P,x));
f2 = simple( - nu*diff(v,x,2) - nu*diff(v,y,2) + w*u + diff(P,y));

f1, f2

%% Convection form
% f1 = simple( - mu*diff(u,x,2) - mu*diff(u,y,2) + u*diff(u,x) + v*diff(u,y) + diff(p,x));
% f2 = simple( - mu*diff(v,x,2) - mu*diff(v,y,2) + u*diff(v,x) + v*diff(v,y) + diff(p,y));



end