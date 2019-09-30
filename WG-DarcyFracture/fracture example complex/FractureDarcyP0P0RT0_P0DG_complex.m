function FractureDarcyP0P0RT0_P0DG_complex(refineMeshStep, cases)

ite = refineMeshStep+1;

%% Mesh
data = readMeshFromComplex(refineMeshStep);
node = data.node; elem = data.elem;
elem2edge = data.elem2edge; edge = data.edge; brother = data.brother;size(brother)
lefBoun = data.lefBoun; rigBoun = data.rigBoun;
botBoun = data.botBoun; topBoun = data.topBoun;
NT = size(elem,1); NE = size(edge,1); NGamma = 103*2^refineMeshStep; % size(brother,1)

fprintf('#Number of elements is: %8.0u\n', NT);
fprintf('#Number of edges is: %8.0u\n', NE);
fprintf('#Number of fractured elements: %8.0u\n', NGamma);

fracLength = sqrt(sum((node(edge(brother(:,1),1),:) - node(edge(brother(:,1),2),:)).^2, 2));

%% Parameter
K = 1e4*ones(NGamma,1);       K(35*2^refineMeshStep+1:64*2^refineMeshStep) = 1e-4; 
lGamma = 1e-4*ones(NGamma,1); psi = 0.75; eitaGamma = lGamma./K;

% System
[A, B, C1] = assembleStiffnessMatrixForPrimalDarcyP0P0RT0(node, elem, elem2edge, NT, NE);
[C2, D1, D2, D3] = assembleStiffnessMatrixInterfaceTermP0P0RT0_complex(NE, NGamma, psi, eitaGamma, brother, fracLength);

% if ite==1
%     DF = TPFA_complex(NGamma, lGamma, K, fracLength, ite);
% elseif ite>1
%     DF = TPFA_complex_2(NGamma, lGamma, K, fracLength, ite);
% end
DF = TPFA_complex_refineMeshStep(NGamma, lGamma, K, fracLength, refineMeshStep);

temp = sparse(NT, NGamma);
leftMatrix = [A      B      temp;
              B'     C1+C2  D1;
              temp'  D2     D3+DF];
rightTerm = zeros(NT+NE+NGamma,1);


%% Solve
uh = zeros(NT+NE+NGamma,1);

if strcmp(cases,'TopToBottom')
    fixed1 = topBoun; fixed2 = botBoun;
    uh(NT+fixed1) = 4; uh(NT+fixed2) = 1;
elseif strcmp(cases,'LeftToRight')
    fixed1 = lefBoun; fixed2 = rigBoun;
    uh(NT+fixed1) = 4; uh(NT+fixed2) = 1;
end

free = setdiff(1:NT+NE+NGamma, NT+[fixed1;fixed2]);
rightTerm = rightTerm - leftMatrix*uh;
% uh(free) = leftMatrix(free,free) \ rightTerm(free); 
uh(free) = amg(leftMatrix(free,free), rightTerm(free));

% [~,p] = chol(leftMatrix(free,free)); 
% 
% matrix0 = leftMatrix(free,free);
% dif = matrix0 - matrix0';
% max(abs(dif(:)))
% 
% norm(rightTerm(free) - leftMatrix(free,free)*uh(free)) / norm(rightTerm(free))
% [condest(leftMatrix(free,free)), nnz(leftMatrix(free,free))/prod(size(leftMatrix(free,free)))]

%% Eliminate interior dof
% invA = A \ speye(NT);
% leftMatrix2 = [C1+C2-B'*invA*B  D1;
%               D2               D3+DF];
% free2 = setdiff(1:NE+NGamma, [fixed1;fixed2]);
% fprintf('#Degrees of freedom are: %8.0u\n', length(free2));
% fprintf('#The sparsity of matrix is: %8.4e\n', nnz(leftMatrix2(free2,free2)) / prod(size(leftMatrix2(free2,free2))));

%% Figure for matrix
ph0 = uh(1:NT);
figure; pdesurf(node', elem', ph0');  shading flat; view(2); colormap('default');
axis equal; axis off; caxis([min(uh) max(uh)]);
hcb = colorbar; set(get(hcb,'Title'),'String','Pressure');

phFrac = uh(NT+NE+1:end); fracXAndPCell = cell(NGamma,1);
for i = 1:NGamma
    bp = node(edge(brother(i,1),1),:); ep = node(edge(brother(i,1),2),:);
    fracXAndPCell{i} = [bp(1) bp(2) phFrac(i);
                        ep(1) ep(2) phFrac(i)];
end

Points = node; T = elem; P = ph0;
if strcmp(cases,'TopToBottom')
    if refineMeshStep==0
        save('dataComplex_Matrix_TopToBottom_1.mat', 'T', 'P', 'Points'); save('dataComplex_Fractu_TopToBottom_1.mat', 'fracXAndPCell');
    elseif refineMeshStep ==1
        save('dataComplex_Matrix_TopToBottom_2.mat', 'T', 'P', 'Points'); save('dataComplex_Fractu_TopToBottom_2.mat', 'fracXAndPCell');
    elseif refineMeshStep==2
        save('dataComplex_Matrix_TopToBottom_3.mat', 'T', 'P', 'Points'); save('dataComplex_Fractu_TopToBottom_3.mat', 'fracXAndPCell');
    end
elseif strcmp(cases,'LeftToRight')
    if refineMeshStep==0
        save('dataComplex_Matrix_LeftToRight_1.mat', 'T', 'P', 'Points'); save('dataComplex_Fractu_LeftToRight_1.mat', 'fracXAndPCell');
    elseif refineMeshStep==1
        save('dataComplex_Matrix_LeftToRight_2.mat', 'T', 'P', 'Points'); save('dataComplex_Fractu_LeftToRight_2.mat', 'fracXAndPCell');
    elseif refineMeshStep==2
        save('dataComplex_Matrix_LeftToRight_3.mat', 'T', 'P', 'Points'); save('dataComplex_Fractu_LeftToRight_3.mat', 'fracXAndPCell');
    end
end
%% figure pressure using Tecplot
% DGM = zeros(3*NT,2); newPh0 = zeros(3*NT,1);
% for i = 1:NT
%     DGM(3*i-2,:) = node(elem(i,1),:);
%     DGM(3*i-1,:) = node(elem(i,2),:);
%     DGM(3*i,:)   = node(elem(i,3),:);
%     newPh0(3*i-2:3*i) = ph0(i);
% end
% DGT = reshape(1:3*NT, 3, NT)';
% 
% figure; trisurf(DGT, DGM(:,1), DGM(:,2), newPh0); shading interp; view(2);
% colorbar;
% f1 = fopen('FractureEquationMatrixPressureComplexLeftToRight.dat','w'); fprintf(f1,'TITLE = ""\n'); fprintf(f1,'VARIABLES = "X" "Y" "U"  \n'); 
% fprintf(f1,'ZONE N=%d,E=%d\n', 3*NT, NT); fprintf(f1,'DATAPACKING=POINT,ZONETYPE=FETRIANGLE\n');
% for i = 1:3*NT
%    fprintf(f1, '%f\t%f\t%f\t\r\n', DGM(i,1), DGM(i,2), newPh0(i));  
% end
% fprintf('\n');
% for n = 1:NT
%    fprintf(f1, '%d\t%d\t%d\r\n', 3*n-2, 3*n-1, 3*n); 
% end
% fclose(f1);
end