
clear clc
syms x y lambda mu n1 n2

% u1 = x.^2*y.^2+exp(-y)
% u2 = -2*x*y.^3/3+2-pi*sin(pi*x)

u1 = cos(pi*x).*cos(pi*y);  
u2 = sin(pi*x).*sin(pi*y);  

u1_x = simplify(diff(u1,x,1))
u1_y = simplify(diff(u1,y,1))
u2_x = simplify(diff(u2,x,1))
u2_y = simplify(diff(u2,y,1))

epsilon = simplify( 0.5*[ u1_x + u1_x,  u1_y + u2_x;
                          u2_x + u1_y,  u2_y + u2_y,])
            
sigma = simplify( lambda*(u1_x+u2_y)*eye(2)+2*mu*epsilon )

f1 = simplify(-(diff(sigma(1,1),x,1)+diff(sigma(1,2),y,1)))
f2 = simplify(-(diff(sigma(2,1),x,1)+diff(sigma(2,2),y,1)))


% n =[n1,n2]';
n =[1,0]';
g_N=simplify( sigma*n)