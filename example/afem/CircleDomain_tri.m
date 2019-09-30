function CircleDomain_tri()

[node,elem] = circlemesh(0,0,1,1/4);
uniformrefine

showmesh(node,elem);
axis on 
save('Lshape2_tri_64','node','elem')

end