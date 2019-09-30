function FractureDarcyP0P0RT0_P0DG_geiger(n, KGammaN)


%% Mesh
[node, elem, elem2edge, edge, NE, NT, freeEdge, leftNeumann, ...
 bottomNeumann, topNeumann, Dirichlet, frac1, frac2, frac3, ...
 frac4, frac5, frac6, N, h] = classifyIntersectionEdge(n);

% figure; trisurf(elem, node(:,1), node(:,2), zeros(size(node,1),1), 'FaceColor','w', 'EdgeColor','k', 'LineWidth', 0.3); 
% view(2); axis equal; axis off; hold on;
% x = cell(6); y = cell(6);
% x{1} = (0:0.01:1)'; y{1} = 0.5*ones(length(x{1}),1);
% x{2} = y{1}; y{2} = x{1};
% x{3} = (0.5:0.01:1)'; y{3} = 0.75*ones(length(x{3}),1);
% x{4} = y{3}; y{4} = x{3};
% x{5} = (0.5:0.01:0.75)'; y{5} = 0.625*ones(length(x{5}),1);
% x{6} = y{5}; y{6} = x{5};
% for i = 1:6
%     plot(x{i}, y{i}, 'LineWidth', 1.8); hold on;
% end
% x = 0: 0.01: 0.9; y = x + 0.1;
% plot(x,y,'r'); hold off; 

%% Parameters
KGammaTau = KGammaN; lGamma = 1e-4; psi = 0.75; eitaGamma = lGamma/KGammaN;


%% System
[A, B, C1, area] = assembleStiffnessMatrixForPrimalDarcyP0P0RT0(node, elem, elem2edge, NT, NE);
[C2, D1, D2, D3, NGamma] = assembleStiffnessMatrixInterfaceTermP0P0RT0_geiger(NE, h, frac1, ...
                                             frac2, frac3, frac4, frac5, frac6, psi, eitaGamma);
[DF, b3] = TPFA_geiger(NGamma, h, n, lGamma, KGammaTau);
% left Neumann
b2 = zeros(NE,1);
for i = 1:length(leftNeumann)
    j = leftNeumann(i); b2(j) = h;
end

temp = sparse(NT, NGamma);
leftMatrix = [A      B      temp;
              B'     C1+C2  D1;
              temp'  D2     D3+DF];
rightTerm = [zeros(NT,1); b2; b3];
clear b2 b3 temp


%% Solve
fixed = NT+Dirichlet; free = setdiff(1:NT+NE+NGamma, fixed);
uh = zeros(NT+NE+NGamma,1); uh(fixed) = 1;

rightTerm = rightTerm - leftMatrix*uh;
% uh(free) = leftMatrix(free,free) \ rightTerm(free);
uh(free) = amg(leftMatrix(free,free), rightTerm(free));

% freeMatrix = leftMatrix(free,free); 
matrix0 = leftMatrix(free,free); [~,p] = chol(matrix0); p
dif = matrix0 - matrix0';
max(abs(dif(:)))

[max(abs(sum(abs(leftMatrix(free,free))))), condest(leftMatrix(free,free)), nnz(leftMatrix(free,free))/prod(size(leftMatrix(free,free)))]

%% Figure for matrix
ph0 = uh(1:NT); phb = uh(NT+1:NT+NE);
figure; pdesurf(node', elem', ph0');  view(2); 
shading flat; view(2); colormap('default');
axis equal; axis off; caxis([min(uh) max(uh)]);
hcb = colorbar; set(get(hcb,'Title'),'String','Pressure');

maxValue = max(ph0); minValue = min(ph0);
step = (maxValue - minValue) / 5;
colorbarValue = minValue: step: maxValue; 
% digits(2); colorbarValue(1) = vpa(colorbarValue(1)); colorbarValue(end) = vpa(colorbarValue(end));
% digits(3); colorbarValue(2:end-1) = vpa(colorbarValue(2:end-1));  colorbarValue
set(hcb,'YTick', colorbarValue);
% t = get(hcb,'YTickLabel');
% set(hcb,'YTickLabel',t);


%% Figure for velocity field
center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
discreteWeakGradient = computeWeakGradientTriangularMesh(center(:,1), center(:,2), node, elem, elem2edge, area, ph0, phb);
figure;
discreteWeakGradient = - discreteWeakGradient;
quiver(center(:,1), center(:,2), discreteWeakGradient(:,1), discreteWeakGradient(:,2));
set(gca, 'fontsize', 8);
axis equal; axis off; view(2);

%% Figure for fracture
phFrac = uh(NT+NE+1:end);
brother = [frac1; frac2; frac3; frac4; frac5; frac6];
fracXAndPCell = cell(NGamma,1);
for i = 1:NGamma
    bp = node(edge(brother(i,1),1),:); ep = node(edge(brother(i,1),2),:);
    fracXAndPCell{i} = [bp(1) bp(2) phFrac(i);
                        ep(1) ep(2) phFrac(i)];
end

Points = node; T = elem; P = ph0;
if KGammaN==1e4
    save('dataGeigerMatrixConductive.mat', 'T', 'P', 'Points'); 
    save('dataGeigerFractureConductive.mat', 'fracXAndPCell');
elseif KGammaN==1e-4
    save('dataGeigerMatrixBlocking.mat', 'T', 'P', 'Points'); 
    save('dataGeigerFractureBlocking.mat', 'fracXAndPCell');
end

center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;

if KGammaN==1e4
    % figure across the longest vertical frature x=0.5
    open geiger_conductive_longest_vertical.fig
    hold on
    phFrac_vertical_fracture = phFrac(N+1:2*N);
    yy = zeros(1,6*N); phFrac_vertical_fractureTemp = zeros(1,6*N);
    phb = uh(NT+1:NT+NE);
    for i = 1:N
        yy((i-1)*6+1:i*6) = (i-1)*h: h/5: i*h; %  + ph0(frac2(i,4)))/2
%         phFrac_vertical_fractureTemp((i-1)*6+1:i*6) = phb(frac2(i,1)) * ones(1,6);
        phFrac_vertical_fractureTemp((i-1)*6+1:i*6) = phFrac_vertical_fracture(i) * ones(1,6);
    end
    plot(yy, phFrac_vertical_fractureTemp, 'Color', [0 0.5 1], 'LineWidth', 1);
    xlabel('arc length', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    ylabel('pressure', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    set(gca, 'FontSize', 12);
    legend('reference', 'Box-DFM', 'CC-DFM', 'EDFM', 'mortar-DFM', 'P-XFEM', 'D-XFEM', 'WG');
    hold off
    clear phFrac_vertical_fracture yy phFrac_vertical_fractureTemp
    
    % figure across y=0.7
    [temp1, unused] = find(center(:,2)>0.6875);
    [temp2, unused] = find(center(:,2)<0.71875);
    temp3 = intersect(temp1, temp2);
    [unused, ix] = sort(center(temp3, 1));
    temp = temp3(ix); % across positive direction of x-axis
    open geiger_conductive_horizontal_line.fig
    hold on
    xx = zeros(1,12*N); ph_horizontal_matrix = ph0(temp); 
    ph_horizontal_matrixTemp = zeros(1,12*N);
    for i = 1:N
        xx((2*i-2)*6+1:(2*i-1)*6) = (i-1)*h: 0.0125/5: 0.0125+(i-1)*h;
        ph_horizontal_matrixTemp((2*i-2)*6+1:(2*i-1)*6) = ph_horizontal_matrix(2*i-1) * ones(1,6);
        xx((2*i-1)*6+1:2*i*6) = 0.0125+(i-1)*h: (h-0.0125)/5: i*h;
        ph_horizontal_matrixTemp((2*i-1)*6+1:2*i*6) = ph_horizontal_matrix(2*i) * ones(1,6);
    end
    plot(xx, ph_horizontal_matrixTemp, 'Color', [0 0.5 1], 'LineWidth', 1); 
    xlabel('arc length', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    ylabel('pressure', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    set(gca, 'FontSize', 12);
    legend('reference', 'Box-DFM', 'CC-DFM', 'EDFM', 'mortar-DFM', 'P-XFEM', 'D-XFEM', 'WG');
    hold off
    clear temp1 temp2 temp3 ix temp xx ph_horizontal_matrix ph_horizontal_matrixTemp
end

if KGammaN==1e-4
    % figure across the line (0,0.1)-(0.9,1)
    left = floor(0.1/h); right = ceil(0.1/h);
    % then, find all elements between y=x+left*h and y=x+right*h
    flag = center(:,2) - center(:,1);
    [temp1, unused] = find(flag>left*h);
    [temp2, unused] = find(flag<right*h);
    temp3 = intersect(temp1, temp2);
    [unused, ix] = sort(center(temp3, 2));
    temp = temp3(ix); 
    open geiger_blocking_diag.fig
    hold on
    hodd = sqrt(2) * abs(4*h-0.1); heven = sqrt(2) * abs(3*h-0.1);
    ss = zeros(1,11*29 + 3*28); ph0ss = zeros(1,11*29 + 3*28);
    tempLength = 0; tempN = 0;
    for i = 1:length(temp)
        if mod(i,2)==1 % i is odd number
            ss(tempN+1:tempN+11) = tempLength: hodd/10: tempLength+hodd;
            ph0ss(tempN+1:tempN+11) = ph0(temp(i)) * ones(1,11);
            tempLength = tempLength + hodd; tempN = tempN + 11;
        elseif mod(i,2)==0 % i is even number
            ss(tempN+1:tempN+3) = tempLength: heven/2: tempLength+heven;
            ph0ss(tempN+1:tempN+3) = ph0(temp(i)) * ones(1,3);
            tempLength = tempLength + heven; tempN = tempN + 3;
        end
    end
    plot(ss, ph0ss, 'Color', [0 0.5 1], 'LineWidth', 1); hold on;
    xlabel('arc length', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    ylabel('pressure', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    set(gca, 'FontSize', 12);
    legend('reference', 'Box-DFM', 'CC-DFM', 'EDFM', 'mortar-DFM', 'P-XFEM', 'D-XFEM', 'WG');
    hold off
    clear left right flag temp1 temp2 temp3 ix temp hodd heven ss ph0ss tempLength
end
fprintf('#Max value is: %8.4e\n', max(uh));
fprintf('#Min value is: %8.4e\n', min(uh));
fprintf('#Mesh size is: %8.0u\n', N);
fprintf('#Number of elements is: %8.0u\n', NT);
fprintf('#Number of edges is: %8.0u\n', NE);

%% figure pressure using Tecplot
% DGM = zeros(3*NT,2); newPh0 = zeros(3*NT,1);
% for i = 1:NT
%     DGM(3*i-2,:) = node(elem(i,1),:);
%     DGM(3*i-1,:) = node(elem(i,2),:);
%     DGM(3*i,:)   = node(elem(i,3),:);
%     newPh0(3*i-2:3*i) = ph0(i);
% end
% 
% f1 = fopen('FractureEquationMatrixPressure.dat','w'); fprintf(f1,'TITLE = ""\n'); fprintf(f1,'VARIABLES = "X" "Y" "U"  \n'); 
% fprintf(f1,'ZONE N=%d,E=%d\n', 3*NT, NT); fprintf(f1,'DATAPACKING=POINT,ZONETYPE=FETRIANGLE\n');
% for i = 1:3*NT
%    fprintf(f1, '%f\t%f\t%f\t\r\n', DGM(i,1), DGM(i,2), newPh0(i));  
% end
% fprintf('\n');
% for n = 1:NT
%    fprintf(f1, '%d\t%d\t%d\r\n', 3*n-2, 3*n-1, 3*n); 
% end
% fclose(f1);

%%
% normest(leftMatrix(free,free),1e-16) * normest(leftMatrix(free,free) \ speye(length(free),length(free)),1e-16)
% InfNormRow = zeros(NT+NE+NGamma,1);
% for i = 1:NT+NE+NGamma
%     InfNormRow(i) = norm(leftMatrix(i,:), inf);
% end
% PrecondMatrix = spdiags(1./InfNormRow, 0, NT+NE+NGamma, NT+NE+NGamma);
% leftMatrix = PrecondMatrix*leftMatrix;
% fprintf('#Condition number of matrix is: %8.4e\n', cond(leftMatrix(free,free)+zeros(length(free),length(free)),2));
% fprintf('#The sparsity of matrix is: %8.4e\n', nnz(leftMatrix(free,free)) / prod(size(leftMatrix(free,free))));

% clear fixed free leftMatrix rightTerm
% save leftMatrix 

% if leftMatrix==leftMatrix'
%     disp('1');
% else
%     disp('0');
% end
% %% Eliminate interior dof
% invA = A \ speye(NT);
% leftMatrix2 = [C1+C2-B'*invA*B  D1;
%                D2               D3+DF];
% free2 = setdiff(1:NE+NGamma, Dirichlet);
% InfNormRow = zeros(NE+NGamma,1);
% for i = 1:NE+NGamma
%     InfNormRow(i) = norm(leftMatrix2(i,:), inf);
% end
% 
% PrecondMatrix = spdiags(1./InfNormRow, 0, NE+NGamma, NE+NGamma);
% leftMatrix3 = PrecondMatrix * leftMatrix2;
% 
% fprintf('#Degrees of freedom are: %8.0u\n', length(free2));
% fprintf('#The sparsity of matrix is: %8.4e\n', nnz(leftMatrix2(free2,free2)) / prod(size(leftMatrix2(free2,free2))));
% fprintf('#Condition number of matrix is: %8.4e\n', cond(leftMatrix2+zeros(NE+NGamma,NE+NGamma),2));

end