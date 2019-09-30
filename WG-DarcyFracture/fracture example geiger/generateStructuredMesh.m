function data = generateStructuredMesh(n)

node1 = [0,0; 0.5,0; 0.5,0.5; 0,0.5]; elem1 = [2,3,1; 4,1,3];
node2 = [0.5,0; 1,0; 1,0.5; 0.5,0.5]; elem2 = [2,3,1; 4,1,3];
node3 = [0,0.5; 0.5,0.5; 0.5,1; 0,1]; elem3 = [2,3,1; 4,1,3];

node4 = [0.75,0.5; 1,0.5; 1,0.75; 0.75,0.75]; elem4 = [2,3,1; 4,1,3];
node5 = [0.75,0.75; 1,0.75; 1,1; 0.75,1]; elem5 = [2,3,1; 4,1,3];
node6 = [0.5,0.75; 0.75,0.75; 0.75,1; 0.5,1]; elem6 = [2,3,1; 4,1,3];

node7 = [0.5,0.5; 0.625,0.5; 0.625,0.625; 0.5,0.625]; elem7 = [2,3,1; 4,1,3];
node8 = [0.625,0.5; 0.75,0.5; 0.75,0.625; 0.625,0.625]; elem8 = [2,3,1; 4,1,3];
node9 = [0.625,0.625; 0.75,0.625; 0.75,0.75; 0.625,0.75]; elem9 = [2,3,1; 4,1,3];
node10 = [0.5,0.625; 0.625,0.625; 0.625,0.75; 0.5,0.75]; elem10 = [2,3,1; 4,1,3];

% n = 0;
for i = 1:n+2
    [node1, elem1] = uniformrefine(node1, elem1);
    [node2, elem2] = uniformrefine(node2, elem2);
    [node3, elem3] = uniformrefine(node3, elem3);
end

for i = 1:n+1
    [node4, elem4] = uniformrefine(node4, elem4);
    [node5, elem5] = uniformrefine(node5, elem5);
    [node6, elem6] = uniformrefine(node6, elem6);
end

for i = 1:n
    [node7, elem7] = uniformrefine(node7, elem7);
    [node8, elem8] = uniformrefine(node8, elem8);
    [node9, elem9] = uniformrefine(node9, elem9);
    [node10, elem10] = uniformrefine(node10, elem10);
end

data = struct('node1',node1, 'elem1',elem1, 'node2',node2, 'elem2',elem2, ...
              'node3',node3, 'elem3',elem3, 'node4',node4, 'elem4',elem4, ...
              'node5',node5, 'elem5',elem5, 'node6',node6, 'elem6',elem6, ...
              'node7',node7, 'elem7',elem7, 'node8',node8, 'elem8',elem8, ...
              'node9',node9, 'elem9',elem9, 'node10',node10, 'elem10',elem10);
          
% figure;
% trisurf(elem1, node1(:,1), node1(:,2), zeros(size(node1,1),1), 'FaceColor', 'w'); view(2); hold on;
% trisurf(elem2, node2(:,1), node2(:,2), zeros(size(node2,1),1)); view(2); hold on;
% trisurf(elem3, node3(:,1), node3(:,2), zeros(size(node3,1),1)); view(2); hold on;
% trisurf(elem4, node4(:,1), node4(:,2), zeros(size(node4,1),1), 'FaceColor', 'w'); view(2); hold on;
% trisurf(elem5, node5(:,1), node5(:,2), zeros(size(node5,1),1)); view(2); hold on;
% trisurf(elem6, node6(:,1), node6(:,2), zeros(size(node6,1),1), 'FaceColor', 'w'); view(2); hold on;
% trisurf(elem7, node7(:,1), node7(:,2), zeros(size(node7,1),1), 'FaceColor', 'r'); view(2); hold on;
% trisurf(elem8, node8(:,1), node8(:,2), zeros(size(node8,1),1), 'FaceColor', 'b'); view(2); hold on;
% trisurf(elem9, node9(:,1), node9(:,2), zeros(size(node9,1),1), 'FaceColor', 'y'); view(2); hold on;
% trisurf(elem10, node10(:,1), node10(:,2), zeros(size(node10,1),1), 'FaceColor', 'g'); view(2); hold off;
% 
% node = [node1; node2; node3; node4; node5; node6; node7; node8; node9; node10];
% nt1 = max(elem1(:)); nt2 = max(elem2(:)); nt3 = max(elem3(:));
% nt4 = max(elem4(:)); nt5 = max(elem5(:)); nt6 = max(elem6(:));
% nt7 = max(elem7(:)); nt8 = max(elem8(:)); nt9 = max(elem9(:));
% elem = [elem1; elem2+nt1; elem3+nt1+nt2; elem4+nt1+nt2+nt3;
%         elem5+nt1+nt2+nt3+nt4;
%         elem6+nt1+nt2+nt3+nt4+nt5;
%         elem7+nt1+nt2+nt3+nt4+nt5+nt6;
%         elem8+nt1+nt2+nt3+nt4+nt5+nt6+nt7;
%         elem9+nt1+nt2+nt3+nt4+nt5+nt6+nt7+nt8;
%         elem10+nt1+nt2+nt3+nt4+nt5+nt6+nt7+nt8+nt9];
% figure; trisurf(elem, node(:,1), node(:,2), zeros(size(node,1),1), 'FaceColor', 'w', 'EdgeColor', 'b'); view(2);
% axis equal; axis off; hold on;
% 
% x = 0:0.05:0.9; plot(x, x+0.1, 'r'); hold off;

end