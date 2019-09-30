function [r1, r2, b] = P0DG_2(n, dis, lGamma, eipsilon, psi, eitaGamma, fGamma, gcrl, gprl, uLef, uRig)

co1 = lGamma*eipsilon; co2 = 4/((2*psi-1)*eitaGamma);

%% 
e = ones(n,1); h = dis/n;
temp = 1/h * spdiags([-e 2*e -e], -1:1, n, n);
temp(1,1) = temp(1,1) + 1/h; temp(n,n) = temp(n,n) + 1/h;
temp = co1 * temp;
r1 = temp; r2 = co2 * h * spdiags(ones(n,1), 0, n, n);

b = zeros(n,1);
for i = 1:n
    beginPoint = [(i-1)*h, 0]; endPoint = [i*h, 0];
    [gcll, gpll] = generate_Gauss_local_1D(gcrl, gprl, beginPoint, endPoint, 0);
    for k = 1:length(gcll)
        b(i,1) = b(i,1) + lGamma * gcll(k) * fGamma(gpll(k));
    end
end
b(1) = b(1) + co1*2*uLef/h;
b(n) = b(n) + co1*2*uRig/h;
end