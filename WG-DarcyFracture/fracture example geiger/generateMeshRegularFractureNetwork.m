function generateMeshRegularFractureNetwork


%% Domain 1
file = fopen('domainNode1.dat'); points = fscanf(file,'%f'); fclose(file);
np = size(points); np = np(1)/3; node1 = zeros(np, 2);
for i = 1:np
    node1(i,1) = points(3*i-1); node1(i,2) = points(3*i);
end
file = fopen('domainElem1.dat'); connecticy = fscanf(file,'%f'); fclose(file);
ne = size(connecticy); ne = ne(1)/9; elem1 = zeros(ne,3);
for i = 1:ne
    elem1(i,1) = connecticy(9*i-6)+1;
    elem1(i,2) = connecticy(9*i-3)+1;
    elem1(i,3) = connecticy(9*i)+1;
end
clear np points connecticy ne 

%% Domain 2
file = fopen('domainNode2.dat'); points = fscanf(file,'%f'); fclose(file);
np = size(points); np = np(1)/3; node2 = zeros(np, 2);
for i = 1:np
    node2(i,1) = points(3*i-1); node2(i,2) = points(3*i);
end
file = fopen('domainElem2.dat'); connecticy = fscanf(file,'%f'); fclose(file);
ne = size(connecticy); ne = ne(1)/9; elem2 = zeros(ne,3);
for i = 1:ne
    elem2(i,1) = connecticy(9*i-6)+1;
    elem2(i,2) = connecticy(9*i-3)+1;
    elem2(i,3) = connecticy(9*i)+1;
end
clear np points connecticy ne 


%% Domain 3
file = fopen('domainNode3.dat'); points = fscanf(file,'%f'); fclose(file);
np = size(points); np = np(1)/3; node3 = zeros(np, 2);
for i = 1:np
    node3(i,1) = points(3*i-1); node3(i,2) = points(3*i);
end
file = fopen('domainElem3.dat'); connecticy = fscanf(file,'%f'); fclose(file);
ne = size(connecticy); ne = ne(1)/9; elem3 = zeros(ne,3);
for i = 1:ne
    elem3(i,1) = connecticy(9*i-6)+1;
    elem3(i,2) = connecticy(9*i-3)+1;
    elem3(i,3) = connecticy(9*i)+1;
end
clear np points connecticy ne 

%% Domain 4
file = fopen('domainNode4.dat'); points = fscanf(file,'%f'); fclose(file);
np = size(points); np = np(1)/3; node4 = zeros(np, 2);
for i = 1:np
    node4(i,1) = points(3*i-1); node4(i,2) = points(3*i);
end
file = fopen('domainElem4.dat'); connecticy = fscanf(file,'%f'); fclose(file);
ne = size(connecticy); ne = ne(1)/9; elem4 = zeros(ne,3);
for i = 1:ne
    elem4(i,1) = connecticy(9*i-6)+1;
    elem4(i,2) = connecticy(9*i-3)+1;
    elem4(i,3) = connecticy(9*i)+1;
end
clear np points connecticy ne 

%% Domain 5
file = fopen('domainNode5.dat'); points = fscanf(file,'%f'); fclose(file);
np = size(points); np = np(1)/3; node5 = zeros(np, 2);
for i = 1:np
    node5(i,1) = points(3*i-1); node5(i,2) = points(3*i);
end
file = fopen('domainElem5.dat'); connecticy = fscanf(file,'%f'); fclose(file);
ne = size(connecticy); ne = ne(1)/9; elem5 = zeros(ne,3);
for i = 1:ne
    elem5(i,1) = connecticy(9*i-6)+1;
    elem5(i,2) = connecticy(9*i-3)+1;
    elem5(i,3) = connecticy(9*i)+1;
end
clear np points connecticy ne 

%% Domain 6
file = fopen('domainNode6.dat'); points = fscanf(file,'%f'); fclose(file);
np = size(points); np = np(1)/3; node6 = zeros(np, 2);
for i = 1:np
    node6(i,1) = points(3*i-1); node6(i,2) = points(3*i);
end
file = fopen('domainElem6.dat'); connecticy = fscanf(file,'%f'); fclose(file);
ne = size(connecticy); ne = ne(1)/9; elem6 = zeros(ne,3);
for i = 1:ne
    elem6(i,1) = connecticy(9*i-6)+1;
    elem6(i,2) = connecticy(9*i-3)+1;
    elem6(i,3) = connecticy(9*i)+1;
end
clear np points connecticy ne 

%% Domain 7
file = fopen('domainNode7.dat'); points = fscanf(file,'%f'); fclose(file);
np = size(points); np = np(1)/3; node7 = zeros(np, 2);
for i = 1:np
    node7(i,1) = points(3*i-1); node7(i,2) = points(3*i);
end
file = fopen('domainElem7.dat'); connecticy = fscanf(file,'%f'); fclose(file);
ne = size(connecticy); ne = ne(1)/9; elem7 = zeros(ne,3);
for i = 1:ne
    elem7(i,1) = connecticy(9*i-6)+1;
    elem7(i,2) = connecticy(9*i-3)+1;
    elem7(i,3) = connecticy(9*i)+1;
end
clear np points connecticy ne 

%% Domain 8
file = fopen('domainNode8.dat'); points = fscanf(file,'%f'); fclose(file);
np = size(points); np = np(1)/3; node8 = zeros(np, 2);
for i = 1:np
    node8(i,1) = points(3*i-1); node8(i,2) = points(3*i);
end
file = fopen('domainElem8.dat'); connecticy = fscanf(file,'%f'); fclose(file);
ne = size(connecticy); ne = ne(1)/9; elem8 = zeros(ne,3);
for i = 1:ne
    elem8(i,1) = connecticy(9*i-6)+1;
    elem8(i,2) = connecticy(9*i-3)+1;
    elem8(i,3) = connecticy(9*i)+1;
end
clear np points connecticy ne 

%% Domain 9
file = fopen('domainNode9.dat'); points = fscanf(file,'%f'); fclose(file);
np = size(points); np = np(1)/3; node9 = zeros(np, 2);
for i = 1:np
    node9(i,1) = points(3*i-1); node9(i,2) = points(3*i);
end
file = fopen('domainElem9.dat'); connecticy = fscanf(file,'%f'); fclose(file);
ne = size(connecticy); ne = ne(1)/9; elem9 = zeros(ne,3);
for i = 1:ne
    elem9(i,1) = connecticy(9*i-6)+1;
    elem9(i,2) = connecticy(9*i-3)+1;
    elem9(i,3) = connecticy(9*i)+1;
end
clear np points connecticy ne 

%% Domain 10
file = fopen('domainNode10.dat'); points = fscanf(file,'%f'); fclose(file);
np = size(points); np = np(1)/3; node10 = zeros(np, 2);
for i = 1:np
    node10(i,1) = points(3*i-1); node10(i,2) = points(3*i);
end
file = fopen('domainElem10.dat'); connecticy = fscanf(file,'%f'); fclose(file);
ne = size(connecticy); ne = ne(1)/9; elem10 = zeros(ne,3);
for i = 1:ne
    elem10(i,1) = connecticy(9*i-6)+1;
    elem10(i,2) = connecticy(9*i-3)+1;
    elem10(i,3) = connecticy(9*i)+1;
end
clear np points connecticy ne 


%% 
% data = struct('node1',node1, 'elem1',elem1, 'node2',node2, 'elem2',elem2, ...
%               'node3',node3, 'elem3',elem3, 'node4',node4, 'elem4',elem4, ...
%               'node5',node5, 'elem5',elem5, 'node6',node6, 'elem6',elem6, ...
%               'node7',node7, 'elem7',elem7, 'node8',node8, 'elem8',elem8, ...
%               'node9',node9, 'elem9',elem9, 'node10',node10, 'elem10',elem10);
% 
% figure;
% trisurf(elem1, node1(:,1), node1(:,2), zeros(size(node1,1),1)); view(2); hold on;
% trisurf(elem2, node2(:,1), node2(:,2), zeros(size(node2,1),1)); view(2); hold on;
% trisurf(elem3, node3(:,1), node3(:,2), zeros(size(node3,1),1)); view(2); hold on;
% trisurf(elem4, node4(:,1), node4(:,2), zeros(size(node4,1),1)); view(2); hold on;
% trisurf(elem5, node5(:,1), node5(:,2), zeros(size(node5,1),1)); view(2); hold on;
% trisurf(elem6, node6(:,1), node6(:,2), zeros(size(node6,1),1)); view(2); hold on;
% trisurf(elem7, node7(:,1), node7(:,2), zeros(size(node7,1),1)); view(2); hold on;
% trisurf(elem8, node8(:,1), node8(:,2), zeros(size(node8,1),1)); view(2); hold on;
% trisurf(elem9, node9(:,1), node9(:,2), zeros(size(node9,1),1)); view(2); hold on;
% trisurf(elem10, node10(:,1), node10(:,2), zeros(size(node10,1),1)); view(2); hold off;

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
figure; trisurf(elem, node(:,1), node(:,2), zeros(size(node,1),1)); view(2);
axis equal; axis off;
end