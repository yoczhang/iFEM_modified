function formAnalyticalSolutions_MHD

syms x y t

%% Example 1
% u1 = y^5 + t^2; u2 = x^5 + t^2; p = 10*(2*x-1)*(2*y-1)*(1+t^2);
% B1 = t^2 + sin(y); B2 = t^2 + sin(x);

%% Example 2
u1 = x^2*(x-1)^2*y*(y-1)*(2*y-1)*cos(t);
u2 = - x*(x-1)*(2*x-1)*y^2*(y-1)^2*cos(t);
p = (2*x-1)*(2*y-1)*cos(t);
B1 = sin(pi*x)*cos(pi*y)*cos(t);
B2 = - cos(pi*x)*sin(pi*y)*cos(t);

%% Example 3
% u1 = exp(t)*cos(y); u2 = 0; B1 = 0; B2 = sin(t)*cos(x); p = - x*cos(y);

%% Form f1, f2, g1, g2
f1 = simple(diff(u1,t) - diff(u1,x,2) - diff(u1,y,2) + diff(p,x) + ...
            u1*diff(u1,x) + u2*diff(u1,y) + (diff(B2,x) - diff(B1,y))*B2);
        
f2 = simple(diff(u2,t) - diff(u2,x,2) - diff(u2,y,2) + diff(p,y) + ...
            u1*diff(u2,x) + u2*diff(u2,y) - (diff(B2,x) - diff(B1,y))*B1);
        
g1 = simple(diff(B1,t) - diff(B1,x,2) - diff(B1,y,2) - diff(u1*B2-u2*B1,y));
g2 = simple(diff(B2,t) - diff(B2,x,2) - diff(B2,y,2) + diff(u1*B2-u2*B1,x));

f1,f2,g1,g2

end