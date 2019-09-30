function [elem2edge, edge, NT, NE, NGamma, fractureLength, brother, leftDiri, rightDiri, leftFractureEdge, rightFractureEdge] = clarifyEdgesReal(node, elem, line)


%% Figure
figure; trisurf(elem, node(:,1), node(:,2), zeros(size(node,1),1), ...
              'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 0.5); view(2);
axis equal; axis off; hold on;
 
for i = 1:size(line,1)
    bp = node(line(i,1),:); ep = node(line(i,2),:);
    fl = (bp(2) + ep(2))/2;
    if fl>0&&fl<600
       x = bp(1): (ep(1) - bp(1))/100: ep(1);
       y = (ep(2) - bp(2)) / (ep(1) - bp(1)) *(x - bp(1)) + bp(2);
       plot(x, y, 'b', 'LineWidth', 1.8); hold on;
    end
end
clear fl bp ep x y
% text(node(line(57,1),1), node(line(57,1),2), int2str(1), 'FontSize',6,'Color','blue'); hold on;
% text(node(line(57,2),1), node(line(57,2),2), int2str(2), 'FontSize',6,'Color','blue');
hold off; 



%%
[elem2edge,edge] = dofedge(elem); 
NE = size(edge,1); NT = size(elem,1); NGamma = size(line,1);

fractureLength = sqrt(sum((node(line(:,1),:) - node(line(:,2),:)).^2, 2));

% modify elem2edge and edge
% brother(i,:) = [number of fractured edge, matrix edge 1 and 2]
brother = zeros(NGamma,3); addEdge = zeros(NGamma,2); 
middleFracture = (node(line(:,1),:) + node(line(:,2),:))/2;
middleAllEdges = (node(edge(:,1),:) + node(edge(:,2),:))/2; 
for i = 1:NGamma
    coor = middleFracture(i,:); coor = [coor(1)*ones(NE,1) coor(2)*ones(NE,1)];
    distance = sqrt(sum((coor - middleAllEdges).^2, 2));
    [interEdge, unused] = find(distance<1e-8);
    
    [ele, flag] = find(elem2edge==interEdge); 

    if ele(1)<ele(2)
        eleModify = ele(2); localIndex = flag(2);
    else
        eleModify = ele(1); localIndex = flag(1);
    end
    addEdge(i,:) = edge(interEdge,:); 
    elem2edge(eleModify, localIndex) = NE + i;
    brother(i,:) = [i, interEdge, NE+i];
end

edge = [edge; addEdge]; NE = NE + NGamma;

middle = (node(edge(:,1),:) + node(edge(:,2),:))/2;
[leftDiri, unused] = find(middle(:,1)==0);
[rightDiri, unused] = find(middle(:,1)==700);

[unused, used] = sort(middleFracture(:,1));
leftFractureEdge = used(1); rightFractureEdge = used(end);

clear addEdge middleFracture middleAllEdges coor distance interEdge 
clear middle ele flag eleModify localIndex
end