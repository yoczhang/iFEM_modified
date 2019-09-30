function test(refine)

node = [0,0; 1,0; 1,1; 0,1]; elem = [1,2,3,4];
for i = 1:refine, [node, elem] = uniformrefinequad(node, elem); end

area = 1/size(elem,1) * ones(size(elem,1),1); 
[node, elem] = generatePerturbingQuadrilateralMesh(node, elem, area, 0.4); 

dim = ones(1,size(elem,1)); elem = mat2cell(elem, dim); 
weakGalerkin_P1P0P0_polygon_test(node, elem);

end