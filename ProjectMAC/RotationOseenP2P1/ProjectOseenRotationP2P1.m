nodec = [0,0; 1,0; 1,1; 0,1]; elemc = [2,3,1; 4,1,3];
bdFlagc = setboundary(nodec,elemc,'Dirichlet'); 
for i = 1:4
    [nodec,elemc,bdFlagc] = uniformrefine(nodec,elemc,bdFlagc);
end
level = 2; nu = 1; data = dataRotation(nu);

%% Form matrices
A1 = cell(level,1); M = cell(level,1); A = cell(level,1); 
B = cell(level,1); Sp = cell(level,1);
ProP = cell(level-1,1); ResP = cell(level-1,1); 
ProU = cell(level-1,1); ResU = cell(level-1,1);
node = cell(level,1); elem = cell(level,1);
node{1} = nodec; elem{1} = elemc;

for k = 1:level-1
    [elem2dofc, edgec, ~] = dofP2(elemc);
    
    [nodef,elemf,bdFlagf,HB] = uniformrefine(nodec,elemc,bdFlagc);
    [elem2doff, edgef, ~] = dofP2(elemf); 
    node{k+1} = nodef; elem{k+1} = elemf;
    
    NTc = size(elemc,1); NTf = size(elemf,1);
    Npc = size(nodec,1); Npf = size(nodef,1); 
    Nuc = Npc + size(edgec,1); Nuf = Npf + size(edgef,1);
    [ProP{k}, ResP{k}] = Pro_Res_HB_way(nodec, nodef, elemf, HB);
    ProU{k} = transferP2red(elemc,elemf); ResU{k} = ProU{k}';
    
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

% We need: A, B, Sp, ProP, ResP, ProU, ResU, f1, f2 

%% Test smoother
% error = 1;
% bigu = [u1; u2; zeros(Npf,1)];
% bigA = [A{level} B{level}'; B{level} sparse(Npf,Npf)];
% bigF = [f1; f2; g];
% while (error>1e-3)
%     bigu = LSC_DGS_rotation_smoother(bigu, A{level}, B{level}, bigF, Sp{level}, nodef, elemf);
%     
%     r = bigF - bigA*bigu;
%     error = norm(r) / norm(bigF)
% end
% 1.5264e-01   3.0238e-01   5.2345e+00
%% Vcycle
% ph = zeros(Npf,1);
% bigA = [A{level} B{level}'; B{level} sparse(Npf,Npf)];
% bigF = [f1; f2; g]; bigu = [u1; u2; ph];
% r = bigF - bigA*bigu;  
% figure, trisurf(elem{level}, node{level}(:,1), node{level}(:,2), r(1:Npf)); shading interp; colorbar; view(3);

ph = zeros(Npf,1);
bigu = cycleRotation(u1, u2, ph,  A, B, f1, f2, g, Sp, ProU, ResU, ProP, ResP, level, node, elem);


%% Direct solver
% bigA = [A{level} B{level}'; B{level} sparse(Npf,Npf)];
% bigF = [f1; f2; g];
% i = 1:2*Nuf+Npf-1;
% bigu = zeros(2*Nuf+Npf,1);
% bigu(i) = bigA(i,i) \ bigF(i);


%% Error
uh = bigu(1:Nuf); vh = bigu(Nuf+1:2*Nuf); ph = bigu(2*Nuf+1:end);
c = sum(mean(ph(elemf),2).*area)/sum(area);
ph = ph - c;

elem2 = [elem2doff(:,[1,6,5]); elem2doff(:,[2,4,6]); elem2doff(:,[3,5,4]); elem2doff(:,[4,5,6])];
node2 = [nodef; (nodef(edgef(:,1),:)+nodef(edgef(:,2),:))/2];
uI = data.exactu(node2); pI = data.exactp(nodef);

[norm(uI(:,1) - uh) norm(uI(:,2) - vh) norm(pI - ph)]
