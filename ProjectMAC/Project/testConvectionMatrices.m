function [D1, D2] = testConvectionMatrices(n)

au = ones(n,n+1); bu = ones(n,n+1); av = ones(n+1,n); bv = ones(n+1,n);

[D1, D2] = formConvectionMatrices(au, bu, av, bv);

dof = n*(n+1);
D1 = D1 + zeros(dof,dof);
D2 = D2 + zeros(dof,dof);

end