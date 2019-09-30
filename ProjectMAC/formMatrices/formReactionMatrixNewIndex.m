function [r1, r2] = formReactionMatrixNewIndex(omegaU, omegaV, n)


dof = n*(n+1);  

%% u
nx = n+1;  ny = n;  LEF = 1:ny;   RIG = ny*(nx-1)+1:ny*nx;
a = false;  a(LEF) = true;  a(RIG) = true; inteNodes = find(~a);
i = mod(inteNodes, n); b = (i==0); i = i + n*b; 
j = ceil(inteNodes/n);
ii = [inteNodes; inteNodes; inteNodes; inteNodes];
k1 = i + (j-2)*(n+1); k2 = i + 1 + (j-2)*(n+1);
k3 = i + (j-1)*(n+1); k4 = i + 1 + (j-1)*(n+1);
jj = [k1; k2; k3; k4];
omegaU = omegaU(:); inteOmegaU = omegaU(inteNodes);
co = [inteOmegaU'; inteOmegaU'; inteOmegaU'; inteOmegaU'];
r1 = sparse(ii(:), jj(:), -co(:)/4, dof, dof); 
clear nx ny LEF RIG a b inteNodes ii k1 k2 k3 k4 jj co

% Modify stencil near boundary
i = 1; j = 2:n; ii = i + n*(j-1); jj1 = i + (j-2)*(n+1); jj2 = i + (j-1)*(n+1);
r1(ii,jj1) = 0; r1(ii,jj2) = 0;
clear i j ii jj1 jj2

i = n; j = 2:n; ii = i + n*(j-1); jj1 = i + 1 + (j-2)*(n+1); jj2 = i + 1 + (j-1)*(n+1);
r1(ii,jj1) = 0; r1(ii,jj2) = 0;
clear i j ii jj1 jj2

%% v
nx = n;  ny = n+1; TOP = 1:ny:ny*(nx-1)+1;  BOT = ny:ny:ny*nx;
a = false;  a(TOP) = true;  a(BOT) = true;  inteNodes = find(~a);
j = floor(inteNodes/(n+1)) + 1; i = mod(inteNodes, n+1);
ii = [inteNodes; inteNodes; inteNodes; inteNodes]; 
k1 = i - 1 + (j-1)*n; k2 = i + (j-1)*n;
k3 = i - 1 + j*n; k4 = i + j*n;
jj = [k1; k2; k3; k4];
omegaV = omegaV(:); inteOmegaV = omegaV(inteNodes);
co = [inteOmegaV'; inteOmegaV'; inteOmegaV'; inteOmegaV'];
r2 = sparse(ii(:), jj(:), co(:)/4, dof, dof);
clear nx ny TOP BOT a inteNodes ii jj k1 k2 k3 k4 co

% Modify stencil near boundary
i = 2:n; j = 1; ii = i + (n+1)*(j-1); jj1 = i - 1 + (j-1)*n; jj2 = i + (j-1)*n;
r2(ii,jj1) = 0; r2(ii,jj2) = 0;
clear i j ii jj1 jj2

i = 2:n; j = n; ii = i + (n+1)*(j-1); jj1 = i - 1 + j*n; jj2 = i + j*n;
r2(ii,jj1) = 0; r2(ii,jj2) = 0;
clear i j ii jj1 jj2
end