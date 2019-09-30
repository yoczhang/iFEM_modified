function [node1,elem1,elem2edge1,edge1,fixedEdge1,freeEdge1,node2,elem2,elem2edge2,edge2,fixedEdge2,freeEdge2,brother,Pb,Tb] = formMeshInformationFracturedDarcy_2(ite, basis_type)

          
%% Omega_1,node1, elem1, NE1, NT1, elem2edge1, edge1 fixedEdge1, freeEdge1
node1 = [0,0; 1,0; 1,1; 0,1]; elem1 = [2,3,1; 4,1,3];
% node1 = [0,0; 0,1/2; 1/2,0; 1/2,1/2; 1,0; 1,1/2]; elem1 = [2,1,4; 3,4,1; 4,3,6; 5,6,3];
for i = 1:ite
    [node1,elem1] = uniformrefine(node1,elem1);
end
[elem2edge1,edge1] = dofedge(elem1);
NE1 = size(edge1,1); NT1 = size(elem1,1); 
s1 = accumarray(elem2edge1(:), 1, [NE1 1]);
fixedEdge1 = find(s1 == 1); freeEdge1 = find(s1 == 2); clear s1

%% Omega_2, node2, elem2, NE2, NT2, elem2edge2, edge2 fixedEdge2, freeEdge2
node2 = [0,-1; 1,-1; 1,0; 0,0]; elem2 = [2,3,1; 4,1,3];
% node2 = [0,-1/2; 0,0; 1/2,-1/2; 1/2,0; 1,-1/2; 1,0]; elem2 = [2,1,4; 3,4,1; 4,3,6; 5,6,3];
for i = 1:ite
    [node2,elem2] = uniformrefine(node2,elem2);
end
[elem2edge2,edge2] = dofedge(elem2);
NE2 = size(edge2,1); NT2 = size(elem2,1); 
s2 = accumarray(elem2edge2(:), 1, [NE2 1]);
fixedEdge2 = find(s2 == 1); freeEdge2 = find(s2 == 2); clear s2

%% Intersect line, m, brother
% divieded by 2^n parts, numbered from left to right
% brother(t,:) are numbers of boundary edges from Omega_1 and Omega_2
middleFixed1 = (node1(edge1(fixedEdge1,1),:) + node1(edge1(fixedEdge1,2),:))/2;
middleFixed2 = (node2(edge2(fixedEdge2,1),:) + node2(edge2(fixedEdge2,2),:))/2;
n = 2^ite; 
% h = 2/n; 
h = 1/n;
brother = zeros(n,2);
for i = 1:n
    cx = h/2 + (i-1)*h; cy = 0; % location of middle of current line
    distance1 = sqrt( (cx*ones(size(middleFixed1,1),1)-middleFixed1(:,1)).^2 + ...
                      (cy*ones(size(middleFixed1,1),1)-middleFixed1(:,2)).^2);
    distance2 = sqrt( (cx*ones(size(middleFixed2,1),1)-middleFixed2(:,1)).^2 + ...
                      (cy*ones(size(middleFixed2,1),1)-middleFixed2(:,2)).^2);
    [temp1,unused] = find(distance1==0);
    [temp2,unused] = find(distance2==0);
    brother(i,:) = [fixedEdge1(temp1) fixedEdge2(temp2)];
end
clear  distance1 distance2 temp1 temp2 

%% FEM mesh, Pb, Tb
if basis_type==1
    Pb = -1 : h : 1; Tb = [1:n; 2:n+1];
elseif basis_type==2
    Pb = -1 : h/2 : 1; Tb = [1:2:2*n-1; 2:2:2*n; 3:2:2*n+1];
elseif basis_type==3
    Pb = -1 : h/3: 1; Tb = [1:3:3*n-2; 2:3:3*n-1; 3:3:3*n; 4:3:3*n+1];
end

end