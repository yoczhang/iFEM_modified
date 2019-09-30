%% Quadrilateral mesh with perturbing interior nodes, eipsilon = 1
node1 = [0,0; 1,0; 1,1; 0,1]; elem1 = [1,2,3,4]; node2 = [0,-1; 1,-1; 1,0; 0,0]; elem2 = [1,2,3,4];
n = 5; eipsilon = 1; error0 = zeros(n,1); errorb = zeros(n,1); H1p = zeros(n,1); H1pGamma = zeros(n,1); maxpGamma = zeros(n,1); dof = zeros(n,1);
for i = 1:n
    [node1, elem1] = uniformrefinequad(node1, elem1); [node2, elem2] = uniformrefinequad(node2, elem2); 
    area1 = 1/size(elem1,1) * ones(size(elem1,1),1); area2 = 1/size(elem2,1) * ones(size(elem2,1),1);
    [node11, elem11] = generatePerturbingQuadrilateralMesh(node1, elem1, area1, 0.4); 
    [node22, elem22] = generatePerturbingQuadrilateralMesh(node2, elem2, area2, 0.4); 

    dim1 = ones(1,size(elem11,1)); elem11 = mat2cell(elem11, dim1);
    dim2 = ones(1,size(elem22,1)); elem22 = mat2cell(elem22, dim2);
    
    [error0(i), errorb(i), H1p(i), H1pGamma(i), maxpGamma(i), dof(i)] = weakGalerkin_P1P0P0_polygon(node11,elem11,node22,elem22, eipsilon);
end
error0Rate = zeros(n,1); errorbRate = zeros(n,1); H1pRate = zeros(n,1); H1pGammaRate = zeros(n,1); maxpGammaRate = zeros(n,1); digits(3);
for i = 2:n
    error0Rate(i) = vpa( log(error0(i-1)/error0(i))/log(2));
    errorbRate(i) = vpa( log(errorb(i-1)/errorb(i))/log(2));
    H1pRate(i) = vpa( log(H1p(i-1)/H1p(i))/log(2));
    H1pGammaRate(i) = vpa( log(H1pGamma(i-1)/H1pGamma(i))/log(2));
    maxpGammaRate(i) = vpa( log(maxpGamma(i-1)/maxpGamma(i))/log(2));
end
display('Table 1: quadrilateral mesh with perturbing interior nodes when eipsilon is 1');
colname = {'#DOF',         'error0', 'rate'  'errorb', 'rate', 'H1p', 'rate',  'H1pGamma', 'rate', 'maxpGamma','rate'};
disptable(colname,dof,[], error0,'%0.4e', error0Rate, [], errorb,'%0.4e', errorbRate, [], ...
            H1p, '%0.4e', H1pRate, [], H1pGamma,'%0.4e', H1pGammaRate, [], maxpGamma,'%0.4e', maxpGammaRate, []);
clear all

%% Quadrilateral mesh with perturbing interior nodes, eipsilon = 1000
node1 = [0,0; 1,0; 1,1; 0,1]; elem1 = [1,2,3,4]; node2 = [0,-1; 1,-1; 1,0; 0,0]; elem2 = [1,2,3,4];
n = 5; eipsilon = 1000; error0 = zeros(n,1); errorb = zeros(n,1); H1p = zeros(n,1); H1pGamma = zeros(n,1); maxpGamma = zeros(n,1); dof = zeros(n,1);
for i = 1:n
    [node1, elem1] = uniformrefinequad(node1, elem1); [node2, elem2] = uniformrefinequad(node2, elem2); 
    area1 = 1/size(elem1,1) * ones(size(elem1,1),1); area2 = 1/size(elem2,1) * ones(size(elem2,1),1);
    [node11, elem11] = generatePerturbingQuadrilateralMesh(node1, elem1, area1, 0.4); 
    [node22, elem22] = generatePerturbingQuadrilateralMesh(node2, elem2, area2, 0.4); 
    dim1 = ones(1,size(elem11,1)); elem11 = mat2cell(elem11, dim1);
    dim2 = ones(1,size(elem22,1)); elem22 = mat2cell(elem22, dim2);
    
    [error0(i), errorb(i), H1p(i), H1pGamma(i), maxpGamma(i), dof(i)] = weakGalerkin_P1P0P0_polygon(node11,elem11,node22,elem22, eipsilon);
end
error0Rate = zeros(n,1); errorbRate = zeros(n,1); H1pRate = zeros(n,1); H1pGammaRate = zeros(n,1); maxpGammaRate = zeros(n,1); digits(3);
for i = 2:n
    error0Rate(i) = vpa( log(error0(i-1)/error0(i))/log(2));
    errorbRate(i) = vpa( log(errorb(i-1)/errorb(i))/log(2));
    H1pRate(i) = vpa( log(H1p(i-1)/H1p(i))/log(2));
    H1pGammaRate(i) = vpa( log(H1pGamma(i-1)/H1pGamma(i))/log(2));
    maxpGammaRate(i) = vpa( log(maxpGamma(i-1)/maxpGamma(i))/log(2));
end
display('Table 2: quadrilateral mesh with perturbing interior nodes when eipsilon is 1000');
colname = {'#DOF',         'error0', 'rate'  'errorb', 'rate', 'H1p', 'rate',  'H1pGamma', 'rate', 'maxpGamma','rate'};
disptable(colname,dof,[], error0,'%0.4e', error0Rate, [], errorb,'%0.4e', errorbRate, [], ...
            H1p, '%0.4e', H1pRate, [], H1pGamma,'%0.4e', H1pGammaRate, [], maxpGamma,'%0.4e', maxpGammaRate, []);
clear all

%% Quadrilateral mesh with perturbing interior nodes, eipsilon = 0.001
node1 = [0,0; 1,0; 1,1; 0,1]; elem1 = [1,2,3,4]; node2 = [0,-1; 1,-1; 1,0; 0,0]; elem2 = [1,2,3,4];
n = 5; eipsilon = 0.001; error0 = zeros(n,1); errorb = zeros(n,1); H1p = zeros(n,1); H1pGamma = zeros(n,1); maxpGamma = zeros(n,1); dof = zeros(n,1);
for i = 1:n
    [node1, elem1] = uniformrefinequad(node1, elem1); [node2, elem2] = uniformrefinequad(node2, elem2); 
    area1 = 1/size(elem1,1) * ones(size(elem1,1),1); area2 = 1/size(elem2,1) * ones(size(elem2,1),1);
    [node11, elem11] = generatePerturbingQuadrilateralMesh(node1, elem1, area1, 0.4); 
    [node22, elem22] = generatePerturbingQuadrilateralMesh(node2, elem2, area2, 0.4); 
    dim1 = ones(1,size(elem11,1)); elem11 = mat2cell(elem11, dim1);
    dim2 = ones(1,size(elem22,1)); elem22 = mat2cell(elem22, dim2);
    
    [error0(i), errorb(i), H1p(i), H1pGamma(i), maxpGamma(i), dof(i)] = weakGalerkin_P1P0P0_polygon(node11,elem11,node22,elem22, eipsilon);
end
error0Rate = zeros(n,1); errorbRate = zeros(n,1); H1pRate = zeros(n,1); H1pGammaRate = zeros(n,1); maxpGammaRate = zeros(n,1); digits(3);
for i = 2:n
    error0Rate(i) = vpa( log(error0(i-1)/error0(i))/log(2));
    errorbRate(i) = vpa( log(errorb(i-1)/errorb(i))/log(2));
    H1pRate(i) = vpa( log(H1p(i-1)/H1p(i))/log(2));
    H1pGammaRate(i) = vpa( log(H1pGamma(i-1)/H1pGamma(i))/log(2));
    maxpGammaRate(i) = vpa( log(maxpGamma(i-1)/maxpGamma(i))/log(2));
end
display('Table 3: quadrilateral mesh with perturbing interior nodes when eipsilon is 0.001');
colname = {'#DOF',         'error0', 'rate'  'errorb', 'rate', 'H1p', 'rate',  'H1pGamma', 'rate', 'maxpGamma','rate'};
disptable(colname,dof,[], error0,'%0.4e', error0Rate, [], errorb,'%0.4e', errorbRate, [], ...
            H1p, '%0.4e', H1pRate, [], H1pGamma,'%0.4e', H1pGammaRate, [], maxpGamma,'%0.4e', maxpGammaRate, []);
clear all