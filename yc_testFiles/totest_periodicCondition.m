% |--- yc test the periodic condition
clear
clc

square = [0,1,0,1];
h = 0.5;


% [node,elem] = squaremesh(square,h);
% showmesh(node,elem); 
% axis on;
% findelem(node,elem);  % plot indices of all triangles
% findnode(node);       % plot indices of all vertices


[node,elem,elem2dof] = squaremeshp(square,h,'x');
% showmesh(node,elem); 
% axis on;
% findelem(node,elem);  % plot indices of all triangles
% findnode(node);   

% set(gcf,'Units','normal'); 
% set(gcf,'Position',[0,0,0.6,0.4]);
% subplot(1,2,1); showmesh(node,elem); findelem(node,elem); findnode(node);
% subplot(1,2,2); showmesh(node,elem); findelem(node,elem); findnodedof(node,elem,elem2dof);



% A1 = sym('a1_%d%d',[3,3]);
% A2 = sym('a2_%d%d',[3,3]);
% A3 = sym('a3_%d%d',[3,3]);
% A4 = sym('a4_%d%d',[3,3]);
% A5 = sym('a5_%d%d',[3,3]);
% A6 = sym('a6_%d%d',[3,3]);
% A7 = sym('a7_%d%d',[3,3]);
% A8 = sym('a8_%d%d',[3,3]);
% 
% AA = {A1, A2, A3, A4, A5, A6, A7, A8};

Nelem = size(elem, 1);
Ndof = max(max(elem));
A_orig = sym(zeros(Ndof,Ndof));
A_peri = sym(zeros(Ndof,Ndof));
F_o = sym(zeros(Ndof,1));
F_p = sym(zeros(Ndof,1));

for el = 1:Nelem
    temp = ['a', int2str(el), '_'];
    temp_f = ['f', int2str(el), '_'];
    for ii = 1:3
        row_o = elem(el, ii);
        row_p = elem2dof(el, ii);
        for jj = 1:3
            col_o = elem(el, jj);
            col_p = elem2dof(el, jj);
            A_orig(row_o, col_o) = A_orig(row_o, col_o) ...
                + sym([temp, int2str(row_o), int2str(col_o)]);
            A_peri(row_p, col_p) = A_peri(row_p, col_p) ...
                + sym([temp, int2str(row_o), int2str(col_o)]);
        end
        F_o(row_o) = F_o(row_o) + sym([temp, int2str(row_o)]);
        F_p(row_p) = F_p(row_p) + sym([temp, int2str(row_o)]);
    end
end


a_o = A_orig;
a_o(1, :) = a_o(1, :) + a_o(7, :);
a_o(2, :) = a_o(2, :) + a_o(8, :);
a_o(3, :) = a_o(3, :) + a_o(9, :);

a_o(:, 1) = a_o(:, 1) + a_o(:, 7);
a_o(:, 2) = a_o(:, 2) + a_o(:, 8);
a_o(:, 3) = a_o(:, 3) + a_o(:, 9);

a_o(7, :) = 0;
a_o(8, :) = 0;
a_o(9, :) = 0;
a_o(:, 7) = 0;
a_o(:, 8) = 0;
a_o(:, 9) = 0;







% |--- end 
disp('end of the file')