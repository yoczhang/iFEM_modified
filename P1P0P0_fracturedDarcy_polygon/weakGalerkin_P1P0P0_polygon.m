function [error0,errorb,H1p,H1pGamma,maxpGamma,dof] = weakGalerkin_P1P0P0_polygon(node1,elem1,node2,elem2, eipsilon)
% date: 10/11/2018
% remark: revision version for the journal, JSC

%% Mesh
[node,elem,elem2edge,edge,flag,area,diameter,centroid,testBoun] = readMesh_WVEM_Stokes_polygon(node1,elem1,node2,elem2);
NT = length(elem); NE = size(edge,1);
edgeLength = sqrt(sum((node(edge(:,1),:) - node(edge(:,2),:)).^2, 2));
psi = 0.75; lGamma = 0.01; eitaGamma = lGamma/eipsilon; rou = 1;

%% Find boundary and interface edges
% interfacesEdges_domainOne, interfacesEdges_domainTwo, otherBounEdges
[bounEdge, unused] = find(testBoun); eps = 1e-10;
edgeCenter = (node(edge(:,1),:) + node(edge(:,2),:)) / 2;
[interfaceEdges, unused] = find(abs(edgeCenter(:,2) - 0) < eps);
isDomainOne = true(length(interfaceEdges), 1); 
isDomainOne(interfaceEdges>NE/2) = false;
interfacesEdges_domainOne = interfaceEdges(isDomainOne); 
interfacesEdges_domainTwo = interfaceEdges(~isDomainOne);
[i1, j1] = sort(edgeCenter(interfacesEdges_domainOne,1)); 
interfacesEdges_domainOne = interfacesEdges_domainOne(j1);
[i2, j2] = sort(edgeCenter(interfacesEdges_domainTwo,1)); 
interfacesEdges_domainTwo = interfacesEdges_domainTwo(j2);
otherBounEdges = setdiff(bounEdge, [interfacesEdges_domainOne; interfacesEdges_domainTwo]);
inteNE = length(interfacesEdges_domainTwo);

%% plot
% elem_temp = cell2mat(elem);
% showmesh(node,elem_temp);
% findelem(node,elem_temp);
% findnode(node);
% findedge(node,edge);

%% Assembling
[A,B,C,Mass,Inte0,rhs] = assembling(NT,NE,node,edge,elem2edge,flag,area,centroid,diameter,edgeLength,rou,@p1,@p2,@f1,@f2);

[I1,I2,M,D,rhs2] = assembleInterfaceMatrix_FVMmatrix(node,edge,lGamma,eipsilon,psi,eitaGamma,...
    NE,inteNE,edgeLength, ...
    interfacesEdges_domainOne,interfacesEdges_domainTwo,@pGamma,@fGamma);
 
bigA = [A                    B     sparse(3*NT,inteNE);
        B'                   C+I1  I2;
        sparse(inteNE,3*NT)  I2'   M+D];
bigF = [rhs; zeros(NE,1); rhs2];

Inteb = computeInterpolation_onFaces(node,edge,NE,@p1,@p2);

%% Solving
dof = 3*NT+NE+inteNE; allDof = (1:dof)'; solu = zeros(dof,1); fixed = 3*NT+otherBounEdges; free = setdiff(allDof, fixed);
solu(fixed) = Inteb(otherBounEdges); bigF = bigF - bigA * solu;
solu(free) = bigA(free,free) \ bigF(free);

ph0 = solu(1:3*NT); phb = solu(3*NT+1:3*NT+NE); phG = solu(3*NT+NE+1:end);

%% Error
diff = [Inte0; Inteb] - [ph0; phb]; Inte_pGamma = pGamma(edgeCenter(interfacesEdges_domainTwo,1));
error0 = sqrt((Inte0-ph0)'*Mass*(Inte0-ph0)); errorb = sqrt(sum((Inteb-phb).^2.*edgeLength)); 
H1p = sqrt(diff'*[A B; B' C]*diff);
H1pGamma = sqrt((phG-Inte_pGamma)'*D*(phG-Inte_pGamma)); maxpGamma = max(abs(phG-Inte_pGamma));

%% Data
    function s = p1(p)
        x = p(:,1); y = p(:,2); s = (1-eipsilon)*exp(x+y) + lGamma*(1-eipsilon)/eipsilon*exp(x) + exp(y) + eitaGamma;
    end
    function s = p2(p)
        x = p(:,1); y = p(:,2); s = (1-eipsilon)*exp(x+y) + exp(y);
    end
    function s = pGamma(x)
        s = (1-eipsilon)*exp(x) + lGamma*(1-eipsilon)/(2*eipsilon)*exp(x) + 1 + eitaGamma/2;
    end
    function s = pGammax(x)
        s = -(exp(x)*(2*eipsilon + lGamma)*(eipsilon - 1))/(2*eipsilon);
    end
    function s = f1(p)
        x = p(:,1); y = p(:,2); s = 2*exp(x + y)*(eipsilon - 1) - exp(y) + (lGamma*exp(x)*(eipsilon - 1))/eipsilon;
    end
    function s = f2(p)
        x = p(:,1); y = p(:,2); s = 2*exp(x + y)*(eipsilon - 1) - exp(y);
    end
    function s = fGamma(x)
        s = (exp(x)*(2*eipsilon + lGamma)*(eipsilon - 1))/2;
    end
end