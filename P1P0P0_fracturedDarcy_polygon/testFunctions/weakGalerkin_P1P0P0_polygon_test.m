function weakGalerkin_P1P0P0_polygon_test(node, elem)

%% Mesh
[elem2edge,edge,flag,area,diameter,centroid,testBoun] = readMesh_WVEM_Stokes_polygon(node,elem);
NT = length(elem); NE = size(edge,1);
edgeLength = sqrt(sum((node(edge(:,1),:) - node(edge(:,2),:)).^2, 2));
rou = 5;

%% Find boundary and interface edges
[bounEdge, unused] = find(testBoun); 
edgeCenter = (node(edge(:,1),:) + node(edge(:,2),:)) / 2;

%% Assembling
[A,B,C,Mass,Inte0,rhs] = assembling(NT,NE,node,edge,elem2edge,flag,area,centroid,diameter,edgeLength,rou,@pre,[],@f,[]);

Inteb = computeInterpolation_onFaces(node,edge,NE,edgeLength,@pre,[]); 
bigA = [A    B;
        B'   C];
bigF = [rhs; zeros(NE,1)];

%% Solving
solu = zeros(3*NT+NE,1); solu(3*NT+bounEdge) = Inteb(bounEdge);
bigF = bigF - bigA * solu; free = setdiff((1:3*NT+NE)', 3*NT+bounEdge);
solu(free) = bigA(free,free) \ bigF(free);
ph0 = solu(1:3*NT); phb = solu(3*NT+1:3*NT+NE);

diff = [Inte0; Inteb] - solu;
[(Inte0-ph0)'*Mass*(Inte0-ph0), sqrt(sum((Inteb-phb).^2.*edgeLength)), sqrt(diff'*bigA*diff)]

%% Data
    function s = pre(p)
        x = p(:,1); y = p(:,2); s = x.*(x-1).*y.*(1-y);
    end
    function s = f(p)
        x = p(:,1); y = p(:,2); s = - 2 * y.*(1-y) - 2*x.*(1-x);
    end
end