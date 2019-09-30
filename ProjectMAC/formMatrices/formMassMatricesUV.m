function [M1, M2] = formMassMatricesUV(n)
% n is mesh size 

lef = (1:n)'; rig = (n^2+1:n*(n+1))';
a = ones(n*(n+1),1); a([lef;rig]) = 0;
M1 = spdiags(a, 0, n*(n+1), n*(n+1));
clear lef rig a

top = (1: n+1: n^2)'; bot = top + n;
a = ones((n+1)*n,1); a([top;bot]) = 0;
M2 = spdiags(a, 0, (n+1)*n, (n+1)*n);
clear top bot a

end