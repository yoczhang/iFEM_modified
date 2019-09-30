function r = formDiffusionMatrixPre(width, n)
% 
h = width / n; dof = n^2;
corner = [1, n, n*(n-1)+1, n^2];
LEF = 2:n-1;   RIG = n*(n-1)+2:n^2-1;
BOT = 2*n:n:(n-1)*n; TOP = n+1:n:(n-2)*n+1;

a = false; boun = [corner, LEF, RIG, BOT, TOP]; 
a(boun) = true; inteNodes = find(~a); 

ii = [inteNodes, inteNodes, inteNodes, inteNodes, ...
      1, 1, n, n, n*(n-1)+1, n*(n-1)+1, n^2, n^2, ...
      LEF, LEF, LEF, ...
      RIG, RIG, RIG, ...
      BOT, BOT, BOT, ...
      TOP, TOP, TOP];
k1 = inteNodes - n; k2 = inteNodes - 1;
k3 = inteNodes + 1; k4 = inteNodes + n;
jj = [k1, k2, k3, k4, ...
      2, 1+n, n-1, 2*n, n*(n-2)+1, n*(n-1)+2, n^2-n, n^2-1, ...
      LEF-1, LEF+1, LEF+n, ...
      RIG-n, RIG-1, RIG+1, ...
      BOT-n, BOT+n, BOT-1, ...
      TOP-n, TOP+n, TOP+1];

ss = -1/h^2*ones(1,length(ii));
r = sparse(ii', jj', ss, dof, dof);
clear ii jj k1 k2 k3 k4 ss

ii = [1, n, n*(n-1)+1, n^2, LEF, RIG, BOT, TOP, inteNodes];
jj = [1, n, n*(n-1)+1, n^2, LEF, RIG, BOT, TOP, inteNodes];
ss = [2/h^2*ones(1,4), 3/h^2*ones(1,4*(n-2)), 4/h^2*ones(1,(n-2)^2)];

r = r + sparse(ii', jj', ss', dof, dof);
clear ii jj ss

end