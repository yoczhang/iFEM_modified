function [node,elem,elem2edge,edge,flag,area,diameter,centroid,testBoun] = readMesh_WVEM_Stokes_polygon(node1,elem1,node2,elem2)

node = [node1; node2]; NT1 = length(elem1); NT2 = length(elem2); NT = NT1 + NT2; elem = cell(NT,1);

%% Get [area, diameter, centroid]
area = zeros(NT,1); diameter = zeros(NT,1); centroid = zeros(NT,2);
totalEdge = []; step = zeros(NT,2); ite = 0;

for n = 1:NT
    if n<=NT1
        vert_ids = elem1{n}; elem{n} = elem1{n};
    else
        vert_ids = elem2{n-NT1}+size(node1,1); elem{n} = elem2{n-NT1}+size(node1,1);
    end
    n_slides = length(vert_ids); verts = node(vert_ids,:);

    step(n,:) = [ite+1, ite+n_slides]; ite = ite + n_slides;
    temp = zeros(n_slides,2); temp(:,1) = vert_ids'; temp(:,2) = vert_ids([2:end,1])';
    totalEdge = [totalEdge; temp];
    
    area_components = verts(:,1).*verts([2:end,1],2) - verts([2:end,1],1).*verts(:,2);
    area(n) = 0.5 * abs(sum(area_components));
    centroid(n,:) = sum((verts + verts([2:end,1],:)) .* repmat(area_components,1,2)) / (6*area(n));
    diam = 0;
    for i = 1:(n_slides-1)
        for j = (i+1):n_slides
            diam = max(diam, norm(verts(i,:) - verts(j,:)));
        end
    end
    diameter(n) = diam;
    clear vert_ids n_slides verts area_components diam temp 
end
clear ite

%% Get [edge, j]
sortedTotalEdge = sort(totalEdge,2); 
flagEdge = sign((totalEdge(:,2) - totalEdge(:,1)) .* (sortedTotalEdge(:,2) - sortedTotalEdge(:,1)));
[edge,i2,j] = myunique(sortedTotalEdge);  % j has the same size as totalEdge
clear sortedTotalEdge totalEdge i2

%% Get [NE, elem2edge]
elem2edge = cell(NT,1); flag = cell(NT,1); NE = size(edge,1); testBoun = zeros(NE,1);

for n = 1:NT
    elem2edge{n} = j(step(n,1):step(n,2))'; flag{n} = flagEdge(step(n,1):step(n,2))';
    testBoun(elem2edge{n}) =  testBoun(elem2edge{n}) + flag{n}';
end

end