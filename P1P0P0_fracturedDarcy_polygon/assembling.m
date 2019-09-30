function [A,B,C,Mass,Inte0,rhs] = assembling(NT,NE,node,edge,elem2edge,flag,area,centroid,diameter,edgeLength,rou,p1,p2,f1,f2)

A = sparse(3*NT,3*NT); B = sparse(3*NT,NE); C = sparse(NE,NE); Mass = sparse(3*NT,3*NT);
Inte0 = zeros(3*NT,1); rhs = zeros(3*NT,1); [lambda, weight] = quadpts(4);

for n = 1:NT 
    localSign = flag{n}'; localArea = area(n); localCen = centroid(n,:); localDia = diameter(n); localEdge = edge(elem2edge{n},:);
    left = node(localEdge(:,1),:); righ = node(localEdge(:,2),:); midd = (left + righ)/2;
    localEdgeLength = edgeLength(elem2edge{n}); ne = length(localEdgeLength); % number of edges in this element
    tau = (righ - left) ./ repmat(localEdgeLength,1,2); normal = [tau(:,2), -tau(:,1)]; % Global tangential and normal vectors
    Normal = repmat(localSign,1,2).*normal; % Actually tangential and normal vectors
    
    M = computeMassMatrix(localEdgeLength, left, righ, localArea, localCen, localDia, Normal);
    dof0 = [3*n-2; 3*n-1; 3*n]; dofb = elem2edge{n}'; 
    
    A00_1 = sum(localEdgeLength); A00_2 = sum((midd(:,1)-localCen(1))./localDia.*localEdgeLength);
    A00_3 = sum((midd(:,2)-localCen(2))./localDia.*localEdgeLength);
    A00_4 = sum((midd(:,1)-localCen(1)).^2./localDia^2.*localEdgeLength); 
    A00_5 = sum((midd(:,1)-localCen(1)).*(midd(:,2)-localCen(2))./localDia^2.*localEdgeLength);
    A00_6 = sum((midd(:,2)-localCen(2)).^2./localDia^2.*localEdgeLength);
    A00 = 1/localDia * [A00_1 A00_2 A00_3; A00_2 A00_4 A00_5; A00_3 A00_5 A00_6];
    
    B0b = - 1/localDia * [localEdgeLength'; (midd(:,1)-localCen(1))'/localDia.*localEdgeLength'; ...
                                            (midd(:,2)-localCen(2))'/localDia.*localEdgeLength'];
    Cbb = 1/localDia * spdiags(localEdgeLength, 0, ne, ne);
    L = 1/localArea * [localEdgeLength'.*Normal(:,1)'; localEdgeLength'.*Normal(:,2)'];
    Cbb = rou * Cbb + localArea * (L') * L; 
    
    trial1 = double(repmat(dof0,1,3));  test1 = double(repmat(dof0',3,1));  A = A + sparse(test1(:), trial1(:), rou*A00(:), 3*NT, 3*NT);
    Mass = Mass + sparse(test1(:), trial1(:), M(:), 3*NT, 3*NT);
    trial2 = double(repmat(dofb',3,1)); test2 = double(repmat(dof0,1,ne)); B = B + sparse(test2(:), trial2(:), rou*B0b(:), 3*NT, NE);
    trial3 = double(repmat(dofb,1,ne)); test3 = double(repmat(dofb',ne,1)); C = C + sparse(test3(:), trial3(:), Cbb(:), NE, NE);
    localInte0 = zeros(3,1); localRhs = zeros(3,1);
    for k = 1:ne
        subLocalArea = abs(det([left(k,:)-localCen; righ(k,:)-localCen])) / 2;
            for j = 1:size(lambda,1)        
                pxy = lambda(j,1)*left(k,:) + lambda(j,2)*righ(k,:) + lambda(j,3)*localCen;
                if localCen(2)>0 % the top subdomain
                    localInte0 = localInte0 + weight(j) * subLocalArea * p1(pxy) * polyBases(pxy,localCen,localDia);
                    localRhs = localRhs + weight(j) * subLocalArea * f1(pxy) * polyBases(pxy,localCen,localDia);
                else
                    localInte0 = localInte0 + weight(j) * subLocalArea * p2(pxy) * polyBases(pxy,localCen,localDia);
                    localRhs = localRhs + weight(j) * subLocalArea * f2(pxy) * polyBases(pxy,localCen,localDia);
                end
            end
    end
    Inte0 = Inte0 + accumarray(dof0, M\localInte0, [3*NT 1]); rhs = rhs + accumarray(dof0, localRhs, [3*NT 1]);
end
    function r = polyBases(s, cen, hE)
        x = s(:,1); y = s(:,2); xE = cen(1); yE = cen(2); r = [1; (x-xE)/hE; (y-yE)/hE];
    end

end