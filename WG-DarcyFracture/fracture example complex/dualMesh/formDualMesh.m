function formDualMesh

data = readMeshFromComplex(0,1);
for i = 1:0, [data.node, data.elem] = uniformrefine(data.node, data.elem); end
[dNode,dElem] = dualmesh(data.node,data.elem);
% hold on;
% figure; 
PolyMshr_PlotMsh(dNode,dElem);

% node1 = [0,0; 1,0; 1,1; 0,1]; elem1 = [2,3,1; 4,1,3];
% node2 = [node1(:,1)+1 node1(:,2)]; elem2 = elem1;
% node3 = [0,1; 1,1; 2,1; 0,2; 1,2; 2,2]; elem3 = [2,5,1; 4,1,5; 3,6,2; 6,5,2];
% for i = 1:3
%     [node1,elem1] = uniformrefine(node1,elem1); [node2,elem2] = uniformrefine(node2,elem2); 
%     [node3,elem3] = uniformrefine(node3,elem3);
% end
% node = [node1; node2; node3]; elem = [elem1; elem2+size(node1,1); elem3+size(node1,1)+size(node2,1)];
% [dNode,dElem] = dualmesh(node,elem);
% % hold on;
% figure; 
% PolyMshr_PlotMsh(dNode,dElem);

end