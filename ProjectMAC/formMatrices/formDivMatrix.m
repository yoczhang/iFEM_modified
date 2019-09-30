function B = formDivMatrix(N, width)

n = N; h = width / n;

nx = n+1;  ny = n;  LEF = 1:ny;   RIG = ny*(nx-1)+1:ny*nx;
a = false;  a(LEF) = true;  a(RIG) = true; inteNodes = find(~a);
i = mod(inteNodes, n); b = (i==0); i = i + n*b; 
j = ceil(inteNodes/n);
ii = [inteNodes inteNodes]'; kk = i + (j-2)*n; jj = [kk kk+n]';
B1 = sparse(ii, jj, [-ones(n*(n-1),1)/h; ones(n*(n-1),1)/h], n*(n+1), n^2);
clear nx ny LEF RIG a inteNodes i j b ii jj kk

nx = n;  ny = n+1; TOP = 1:ny:ny*(nx-1)+1;  BOT = ny:ny:ny*nx;
a = false;  a(TOP) = true;  a(BOT) = true;  inteNodes = find(~a);
j = floor(inteNodes/(n+1)) + 1; i = mod(inteNodes, n+1);
ii = [inteNodes inteNodes]'; kk = n*(j-1) + i; jj = [kk-1 kk]';
B2 = sparse(ii, jj, [ones(n*(n-1),1)/h; -ones(n*(n-1),1)/h], n*(n+1), n^2);
clear nx ny TOP BOT a inteNodes i j ii jj kk

B = [B1; B2];
B = B'; 
clear B1 B2