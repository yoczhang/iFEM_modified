function [A, loadVector] = Poisson1D_2(P, T, basis_type, n, gcrl, gprl, lGamma, eipsilon, psi, eitaGamma, f)

gloDof = length(P); lodof = basis_type + 1; h = 2/n;

%% Stiffness matrix and load vector
A = sparse(gloDof, gloDof); loadVector = zeros(gloDof, 1);

for k = 1:n
    beginPoint = [-1+(k-1)*h, 0]; endPoint = [-1+k*h, 0];
    [gcll, gpll] = generate_Gauss_local_1D(gcrl, gprl, beginPoint, endPoint, 0);
    
    sA = zeros(lodof, lodof); sLoadVector = zeros(lodof, 1);
    for i = 1:length(gcll)
        s1 = localBasis1D(gpll(i), 1, basis_type, beginPoint, endPoint);
        s2 = localBasis1D(gpll(i), 0, basis_type, beginPoint, endPoint);
        s3 = feval(f, gpll(i))*localBasis1D(gpll(i), 0, basis_type, beginPoint, endPoint);
        sA = sA + gcll(i)*(lGamma*eipsilon*s1'*s1 + 4/((2*psi-1)*eitaGamma)*s2'*s2);
        sLoadVector = sLoadVector + gcll(i)*lGamma*s3';
    end
    
    trial = repmat(T(:,k), lodof, 1);
    test = repmat(T(:,k)', lodof, 1);
    test = test(:);
    
    temp = sparse(test, trial, sA(:), gloDof, gloDof);
    A = A + temp;
    
    temp2 = accumarray(T(:,k), sLoadVector, [gloDof 1]);
    loadVector = loadVector + temp2;
end

   
end