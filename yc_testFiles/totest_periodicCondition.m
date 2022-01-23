% |--- yc test the periodic condition
clear

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


for el = 1:Nelem
    temp = ['a', int2str(el), '_'];
    for ii = 1:3
        row_o = elem(el, ii);
        row = elem2dof(el, ii);
        for jj = 1:3
            col_o = elem(el, jj);
            col = elem2dof(el, jj);
            A_orig(row_o, col_o) = A_orig(row_o, col_o) ...
                + sym([temp, int2str(row_o), int2str(col_o)]);
            A_peri(row, col) = A_peri(row, col) ...
                + sym([temp, int2str(row_o), int2str(col_o)]);
        end
    end
end










% |--- end 
disp('end of the file')