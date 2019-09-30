function FractureDarcyP0P0RT0_P0DG_real

%% Mesh
[node, elem, line] = readMeshRealFromBoxDfm;
[elem2edge,edge] = dofedge(elem); 
center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;

%% y = 500, to get xcoor and eleY
Horizontal = node(edge(:,1),2) - node(edge(:,2),2); [flagHorizontal, unused] = find(Horizontal==0);
others = setdiff(1:size(edge,1), flagHorizontal); 

flag = (500 - (node(edge(others,1),2) + node(edge(others,2),2))/2) ./ ((node(edge(others,1),2) - node(edge(others,2),2))/2);
[used1, unused] = find(flag>=-1); [used2, unused] = find(flag<=1);
used = intersect(used1, used2); Used = others(used); 
xcoor = (node(edge(Used,1),1)-node(edge(Used,2),1))/2 .* flag(used) + (node(edge(Used,2),1)+node(edge(Used,1),1))/2;
xcoor = sort(xcoor);
clear Horizontal flagHorizontal others flag used1 used2 used Used unused

yValue = [node(elem(:,1),2) node(elem(:,2),2) node(elem(:,3),2)]';
yValueMax = max(yValue); yValueMax = yValueMax';
yValueMin = min(yValue); yValueMin = yValueMin';
flag_1 = 500 - yValueMax; flag_2 = 500 - yValueMin;
flag = flag_1.*flag_2;
[ele, unused] = find(flag<0); 
[unused, ix] = sort(center(ele,1)); eleY = ele(ix);
clear yValue yValueMax yValueMin flag_1 flag_2 flag ele ix unused

%% x = 625, to get ycoor and eleX
vertical = node(edge(:,1),1) - node(edge(:,2),1); [flagVertical, unused] = find(vertical==0);
others = setdiff(1:size(edge,1), flagVertical); 

flag = (625 - (node(edge(others,1),1) + node(edge(others,2),1))/2) ./ ((node(edge(others,1),1) - node(edge(others,2),1))/2);
[used1, unused] = find(flag>=-1); [used2, unused] = find(flag<=1);
used = intersect(used1, used2); Used = others(used); 
ycoor = (node(edge(Used,1),2)-node(edge(Used,2),2))/2 .* flag(used) + (node(edge(Used,2),2)+node(edge(Used,1),2))/2;
ycoor = sort(ycoor);
clear vertical flagVertical others flag used1 used2 used Used 

xValue = [node(elem(:,1),1) node(elem(:,2),1) node(elem(:,3),1)]';
xValueMax = max(xValue); xValueMax = xValueMax';
xValueMin = min(xValue); xValueMin = xValueMin';
flag_1 = 625 - xValueMax; flag_2 = 625 - xValueMin;
flag = flag_1.*flag_2;
[ele, unused] = find(flag<0); 
[unused, ix] = sort(center(ele,2)); eleX = ele(ix); 
clear xValue xValueMax xValueMin flag_1 flag_2 flag ele unused ix
clear elem2edge edge

%%
[elem2edge, edge, NT, NE, NGamma, fracLength, brother, leftDiri, rightDiri, ...
    leftFractureEdge, rightFractureEdge] = clarifyEdgesReal(node, elem, line);

fprintf('#Number of elements is: %8.0u\n', NT);
fprintf('#Number of edges is: %8.0u\n', NE);
fprintf('#Number of fractured elements: %8.0u\n', NGamma);

%% Parameters
Km = 1e-14; Kf = 1e-8; lGamma = 1e-2; psi = 0.75; eitaGamma = lGamma/Kf;


%% System
[A, B, C1] = assembleStiffnessMatrixForPrimalDarcyP0P0RT0(node, elem, elem2edge, NT, NE);

[C2, D1, D2, D3] = assembleStiffnessMatrixInterfaceTermP0P0RT0_real(NE, NGamma, psi, eitaGamma, brother, fracLength);

[DF, b] = TPFA_real(NGamma, lGamma, Kf, fracLength, line, node, leftFractureEdge, rightFractureEdge);

temp = sparse(NT, NGamma);
leftMatrix = [Km*A    Km*B      temp;
              Km*B'   Km*C1+C2  D1;
              temp'   D2        D3+DF];
rightTerm = [zeros(NT+NE,1); b];


%% Solve 
uh = zeros(NT+NE+NGamma,1);
uh(NT+leftDiri) = 1013250; uh(NT+rightDiri) = 0;

free = setdiff(1:NT+NE+NGamma, NT+[leftDiri; rightDiri]);
rightTerm = rightTerm - leftMatrix*uh;
% uh(free) = leftMatrix(free,free) \ rightTerm(free); 
uh(free) = amg(leftMatrix(free,free), rightTerm(free)); 
matrix0 = leftMatrix(free,free);
[~,p] = chol(matrix0); p
dif = matrix0 - matrix0';
max(abs(dif(:)))

[condest(leftMatrix(free,free)), nnz(leftMatrix(free,free))/prod(size(leftMatrix(free,free))), ...
    NT, NE-length([leftDiri; rightDiri]), NGamma]
%% Figure 
ph0 = uh(1:NT); 
figure; pdesurf(node', elem', ph0'); shading flat; view(2); colormap('default');
axis equal; axis off; caxis([min(uh) max(uh)]);
hcb = colorbar; set(get(hcb,'Title'),'String','Pressure');

% y = 500, using xcoor and eleY
open realYEqual500.fig; hold on;  
xdis = zeros(1,6*length(eleY)); value = zeros(1,6*length(eleY));
for i = 1:length(eleY)
    xdis((i-1)*6+1:i*6) = xcoor(i): (xcoor(i+1)-xcoor(i))/5: xcoor(i+1);
    value((i-1)*6+1:i*6) = ph0(eleY(i)) * ones(1,6);
end
plot(xdis, value, 'Color', [0 0.5 1], 'LineWidth', 1);
xlabel('arc length [m]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
ylabel('pressure [Pa]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'FontSize', 12);
legend('Box-DFM', 'CC-DFM', 'EDFM', 'mortar-DFM', 'Location', 'SouthWest', 'WG');
hold off;
clear xdis value 

% x = 625, using ycoor and eleX
open realXEqual625.fig; hold on; 
ydis = zeros(1,6*length(eleX)); value = zeros(1,6*length(eleX));
for i = 1:length(eleX)
    ydis((i-1)*6+1:i*6) = ycoor(i): (ycoor(i+1)-ycoor(i))/5: ycoor(i+1);
    value((i-1)*6+1:i*6) = ph0(eleX(i)) * ones(1,6);
end
plot(ydis, value, 'Color', [0 0.5 1], 'LineWidth', 1);
xlabel('arc length [m]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
ylabel('pressure [Pa]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'FontSize', 12);
legend('Box-DFM', 'CC-DFM', 'EDFM', 'mortar-DFM', 'Location', 'SouthWest', 'WG');
hold off;
clear ydis value

% plot the locations of y = 500 and x = 625
% open domainFigure.fig
% hold on ; x = 0: 1: 700; y = 500*ones(1,length(x)); plot(x, y, 'b'); 
% y2 = 0: 1: 600; x2 = 625*ones(1,length(y2)); plot(x2, y2, 'b'); hold off;


end