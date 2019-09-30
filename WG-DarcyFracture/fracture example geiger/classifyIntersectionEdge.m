function [node, elem, elem2edge, edge, NE, NT, freeEdge, leftNeumann, bottomNeumann, topNeumann, Dirichlet, ...
         frac1, frac2, frac3, frac4, frac5, frac6, N, h] = classifyIntersectionEdge(n)

%% Load mesh
data = generateStructuredMesh(n);
node1 = data.node1; elem1 = data.elem1;   node2 = data.node2;  elem2 = data.elem2;
node3 = data.node3; elem3 = data.elem3;   node4 = data.node4;  elem4 = data.elem4;
node5 = data.node5; elem5 = data.elem5;   node6 = data.node6;  elem6 = data.elem6;
node7 = data.node7; elem7 = data.elem7;   node8 = data.node8;  elem8 = data.elem8;
node9 = data.node9; elem9 = data.elem9;  node10 = data.node10; elem10 = data.elem10;

node = [node1; node2; node3; node4; node5; node6; node7; node8; node9; node10];
nt1 = max(elem1(:)); nt2 = max(elem2(:)); nt3 = max(elem3(:));
nt4 = max(elem4(:)); nt5 = max(elem5(:)); nt6 = max(elem6(:));
nt7 = max(elem7(:)); nt8 = max(elem8(:)); nt9 = max(elem9(:));
elem = [elem1; elem2+nt1; elem3+nt1+nt2; elem4+nt1+nt2+nt3;
        elem5+nt1+nt2+nt3+nt4;
        elem6+nt1+nt2+nt3+nt4+nt5;
        elem7+nt1+nt2+nt3+nt4+nt5+nt6;
        elem8+nt1+nt2+nt3+nt4+nt5+nt6+nt7;
        elem9+nt1+nt2+nt3+nt4+nt5+nt6+nt7+nt8;
        elem10+nt1+nt2+nt3+nt4+nt5+nt6+nt7+nt8+nt9];
% [node,elem] = uniformrefine(node,elem);
% figure; trisurf(elem, node(:,1), node(:,2), zeros(size(node,1),1)); view(2); axis equal; axis off;
clear node1 node2 node3 node4 node5 node6 node7 node8 node9 node10
clear elem1 elem2 elem3 elem4 elem5 elem6 elem7 elem8 elem9 elem10
clear nt1 nt2 nt3 nt4 nt5 nt6 nt7 nt8 nt9 nt10 data

%% Form elem2edge, edge, NT, NE, fixedEdge, freeEdge for every part
[elem2edge,edge] = dofedge(elem); NE = size(edge,1); NT = size(elem,1); 
s1 = accumarray(elem2edge(:), 1, [NE 1]);
fixedEdge = find(s1 == 1); freeEdge = find(s1 == 2); clear s1


%% 
% coordinates of middle point for every edge
middleFixed = (node(edge(fixedEdge,1),:) + node(edge(fixedEdge,2),:))/2; 

% boundary edges, including Neumann and Dirichlet types
[temp1, unused] = find(middleFixed(:,1)==0);
leftNeumann = fixedEdge(temp1);

[temp2, unused] = find(middleFixed(:,2)==0);
bottomNeumann = fixedEdge(temp2);

[temp3, unused] = find(middleFixed(:,2)==1);
topNeumann = fixedEdge(temp3);

[temp4, unused] = find(middleFixed(:,1)==1);
Dirichlet = fixedEdge(temp4);
clear temp1 temp2 temp3 temp4 unused 

otherEdge = setdiff(fixedEdge, [leftNeumann; bottomNeumann; topNeumann; Dirichlet]);
middleOtherEdge = (node(edge(otherEdge,1),:) + node(edge(otherEdge,2),:))/2; 
middle = (node(edge(:,1),:) + node(edge(:,2),:))/2; 

% partitioned into 2^(n+3) elements in x- and y- directions
N = 2^(n+3); h = 1 / N;

% 1st fracture, frac1 = [bot_edge bot_elem top_edge top_elem], left to
% right
frac1 = zeros(N,4);
[temp, unused] = find(middleOtherEdge(:,2)==0.5);
temp2 = otherEdge(temp);
for i = 1:N
    x = h/2 + (i-1)*h;
    dis = abs(ones(2*N,1)*x - middle(temp2,1));
    [used, unused] = find(dis==0); j1 = temp2(used(1)); j2 = temp2(used(2));
    [ele1, unused] = find(elem2edge==j1);
    [ele2, unused] = find(elem2edge==j2); 
    center1 = (node(elem(ele1,1),:) + node(elem(ele1,2),:) + node(elem(ele1,3),:))/3;
    center2 = (node(elem(ele2,1),:) + node(elem(ele2,2),:) + node(elem(ele2,3),:))/3;
    if center1(2)>center2(2)
        frac1(i,:) = [j2 ele2 j1 ele1];
    else
        frac1(i,:) = [j1 ele1 j2 ele2];
    end
    clear x dis used j1 j2 ele1 ele2 center1 center2
end
clear temp unused temp2

% 2nd fracture, frac2 = [left_edge left_elem right_edge right_elem], bot to
% top
frac2 = zeros(N,4);
[temp, unused] = find(middleOtherEdge(:,1)==0.5);
temp2 = otherEdge(temp);
for i = 1:N
    y = h/2 + (i-1)*h;
    dis = abs(ones(2*N,1)*y - middle(temp2,2));
    [used, unused] = find(dis==0); j1 = temp2(used(1)); j2 = temp2(used(2));
    [ele1, unused] = find(elem2edge==j1);
    [ele2, unused] = find(elem2edge==j2); 
    center1 = (node(elem(ele1,1),:) + node(elem(ele1,2),:) + node(elem(ele1,3),:))/3;
    center2 = (node(elem(ele2,1),:) + node(elem(ele2,2),:) + node(elem(ele2,3),:))/3;
    if center1(1)>center2(1)
        frac2(i,:) = [j2 ele2 j1 ele1];
    else
        frac2(i,:) = [j1 ele1 j2 ele2];
    end
    clear x dis used j1 j2 ele1 ele2 center1 center2
end
clear temp unused temp2    

% 3rd fracture, frac3 = [bot_edge bot_elem top_edge top_elem], left to
% right
frac3 = zeros(N/2,4);
[temp, unused] = find(middleOtherEdge(:,2)==0.75);
temp2 = otherEdge(temp);
for i = 1:N/2
    x = 0.5 + h/2 + (i-1)*h;
    dis = abs(ones(N,1)*x - middle(temp2,1));
    [used, unused] = find(dis==0); j1 = temp2(used(1)); j2 = temp2(used(2));
    [ele1, unused] = find(elem2edge==j1);
    [ele2, unused] = find(elem2edge==j2); 
    center1 = (node(elem(ele1,1),:) + node(elem(ele1,2),:) + node(elem(ele1,3),:))/3;
    center2 = (node(elem(ele2,1),:) + node(elem(ele2,2),:) + node(elem(ele2,3),:))/3;
    if center1(2)>center2(2)
        frac3(i,:) = [j2 ele2 j1 ele1];
    else
        frac3(i,:) = [j1 ele1 j2 ele2];
    end
    clear x dis used j1 j2 ele1 ele2 center1 center2
end
clear temp unused temp2

% 4th fracture, frac2 = [left_edge left_elem right_edge right_elem], bot to
% top
frac4 = zeros(N/2,4);
[temp, unused] = find(middleOtherEdge(:,1)==0.75);
temp2 = otherEdge(temp);
for i = 1:N/2
    y = 0.5 + h/2 + (i-1)*h;
    dis = abs(ones(N,1)*y - middle(temp2,2));
    [used, unused] = find(dis==0); j1 = temp2(used(1)); j2 = temp2(used(2));
    [ele1, unused] = find(elem2edge==j1);
    [ele2, unused] = find(elem2edge==j2); 
    center1 = (node(elem(ele1,1),:) + node(elem(ele1,2),:) + node(elem(ele1,3),:))/3;
    center2 = (node(elem(ele2,1),:) + node(elem(ele2,2),:) + node(elem(ele2,3),:))/3;
    if center1(1)>center2(1)
        frac4(i,:) = [j2 ele2 j1 ele1];
    else
        frac4(i,:) = [j1 ele1 j2 ele2];
    end
    clear x dis used j1 j2 ele1 ele2 center1 center2
end
clear temp unused temp2 

% 5th fracture, frac3 = [bot_edge bot_elem top_edge top_elem], left to
% right
frac5 = zeros(N/4,4);
[temp, unused] = find(middleOtherEdge(:,2)==0.625);
temp2 = otherEdge(temp);
for i = 1:N/4
    x = 0.5 + h/2 + (i-1)*h;
    dis = abs(ones(N/2,1)*x - middle(temp2,1));
    [used, unused] = find(dis==0); j1 = temp2(used(1)); j2 = temp2(used(2));
    [ele1, unused] = find(elem2edge==j1);
    [ele2, unused] = find(elem2edge==j2); 
    center1 = (node(elem(ele1,1),:) + node(elem(ele1,2),:) + node(elem(ele1,3),:))/3;
    center2 = (node(elem(ele2,1),:) + node(elem(ele2,2),:) + node(elem(ele2,3),:))/3;
    if center1(2)>center2(2)
        frac5(i,:) = [j2 ele2 j1 ele1];
    else
        frac5(i,:) = [j1 ele1 j2 ele2];
    end
    clear x dis used j1 j2 ele1 ele2 center1 center2
end
clear temp unused temp2

% 6th fracture, frac2 = [left_edge left_elem right_edge right_elem], bot to
% top
frac6 = zeros(N/4,4);
[temp, unused] = find(middleOtherEdge(:,1)==0.625);
temp2 = otherEdge(temp);
for i = 1:N/4
    y = 0.5 + h/2 + (i-1)*h;
    dis = abs(ones(N/2,1)*y - middle(temp2,2));
    [used, unused] = find(dis==0); j1 = temp2(used(1)); j2 = temp2(used(2));
    [ele1, unused] = find(elem2edge==j1);
    [ele2, unused] = find(elem2edge==j2); 
    center1 = (node(elem(ele1,1),:) + node(elem(ele1,2),:) + node(elem(ele1,3),:))/3;
    center2 = (node(elem(ele2,1),:) + node(elem(ele2,2),:) + node(elem(ele2,3),:))/3;
    if center1(1)>center2(1)
        frac6(i,:) = [j2 ele2 j1 ele1];
    else
        frac6(i,:) = [j1 ele1 j2 ele2];
    end
    clear x dis used j1 j2 ele1 ele2 center1 center2
end
clear temp unused temp2


end