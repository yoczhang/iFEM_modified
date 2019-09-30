nodec = [0,0; 1,0; 1,1; 0,1]; elemc = [2,3,1; 4,1,3];
bdFlagc = setboundary(nodec,elemc,'Dirichlet'); 
level = 7; nu = 1e-3; data = dataRotation2(nu);

%% Form matrices
A1 = cell(level,1); M = cell(level,1); A = cell(level,1); 
B = cell(level,1); Sp = cell(level,1);
ProP = cell(level-1,1); ResP = cell(level-1,1); 
ProU = cell(level-1,1); ResU = cell(level-1,1);

for k = 1:level-1
    [elem2dofc, edgec, ~] = dofP2(elemc);
    
    [nodef,elemf,bdFlagf,HB] = uniformrefine(nodec,elemc,bdFlagc);
    [elem2doff, edgef, ~] = dofP2(elemf); 
    
    NTc = size(elemc,1); NTf = size(elemf,1);
    Npc = size(nodec,1); Npf = size(nodef,1); 
    Nuc = Npc + size(edgec,1); Nuf = Npf + size(edgef,1);
    [ProU{k}, ResU{k}, ProP{k}, ResP{k}] = formProRes(Npc, Npf, HB, Nuc, Nuf, edgef);
    
    [Dlambda,area] = gradbasis(nodec,elemc);
    [A1{k}, M{k}] = formDiffusionMassMatrices(NTc, elem2dofc, Dlambda, area, nodec, elemc, data.diffCo, data.omega, Nuc);
    A{k} = [A1{k} -M{k}; M{k} A1{k}];
    B{k} = divergenceOperator(Nuc, Npc, Dlambda, area, elemc, elem2dofc);
    Sp{k} = formSp(nodec, elemc);
    
    [A{k}, B{k}] = modifyMatrix(Nuc, Npc, bdFlagc, elem2dofc, edgec, A{k}, B{k});
    
    nodec = nodef; elemc = elemf; bdFlagc = bdFlagf;
end
clear nodec elemc edgec elem2dofc bdFlagc NTc Npc Nuc Dlambda area HB

[Dlambda,area] = gradbasis(nodef,elemf);
[A1{level}, M{level}] = formDiffusionMassMatrices(NTf, elem2doff, Dlambda, area, nodef, elemf, data.diffCo, data.omega, Nuf);
A{level} = [A1{level} -M{level}; M{level} A1{level}];

B{level} = divergenceOperator(Nuf, Npf, Dlambda, area, elemf, elem2doff);
Sp{level} = formSp(nodef, elemf);

[f1, f2] = formLoadVector(Nuf, NTf, elem2doff, area, nodef, elemf, data.f);
[u1, u2, f1, f2, g] = modifyDirichletBoundaryCondition(Nuf, Npf, nodef, elem2doff, bdFlagf, edgef, A{level}, B{level}, f1, f2, data.g_D);

[A{level}, B{level}] = modifyMatrix(Nuf, Npf, bdFlagf, elem2doff, edgef, A{level}, B{level});
clear bdFlagf 
clear NTf Dlambda  A1 M 

%%
bigA = [A{level} B{level}'; B{level} sparse(Npf,Npf)];
bigF = [f1; f2; g];
i = 1:2*Nuf+Npf-1;
bigu = zeros(2*Nuf+Npf,1);
bigu(i) = bigA(i,i) \ bigF(i);

uh = bigu(1:Nuf); vh = bigu(Nuf+1:2*Nuf); ph = bigu(2*Nuf+1:end);
c = sum(mean(ph(elemf),2).*area)/sum(area);
ph = ph - c;

elem2 = [elem2doff(:,[1,6,5]); elem2doff(:,[2,4,6]); elem2doff(:,[3,5,4]); elem2doff(:,[4,5,6])];
node2 = [nodef; (nodef(edgef(:,1),:)+nodef(edgef(:,2),:))/2];

vel = sqrt(uh.^2 + vh.^2);
figure, trisurf(elem2, node2(:,1), node2(:,2), vel); shading interp; view(2); colorbar;
figure, trisurf(elemf, nodef(:,1), nodef(:,2), ph); shading interp; view(2); colorbar;

% velocity
z=fopen('NSRotationVelocity.dat','w'); fprintf(z,'TITLE = ""\n'); fprintf(z,'VARIABLES = "X" "Y" "U" "U1" "U2" \n'); 
fprintf(z,'ZONE N=%d,E=%d\n',size(node2,1),size(elem2,1)); fprintf(z,'DATAPACKING=POINT,ZONETYPE=FETRIANGLE\n');
for i = 1:size(node2,1)
   fprintf(z,'%f\t%f\t%f\t%f\t%f\r\n',node2(i,1),node2(i,2),vel(i),uh(i),vh(i));  
end
fprintf('\n');
for n = 1:size(elem2,1)
   fprintf(z,'%d\t%d\t%d\r\n',elem2(n,1),elem2(n,2),elem2(n,3)); 
end
fclose(z);