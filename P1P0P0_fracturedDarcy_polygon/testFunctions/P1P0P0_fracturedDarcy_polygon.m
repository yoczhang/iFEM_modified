function P1P0P0_fracturedDarcy_polygon(refine)

%% Quadrilateral mesh with perturbing interior nodes
node1 = [0,0; 1,0; 1,1; 0,1]; elem1 = [1,2,3,4];
node2 = [0,-1; 1,-1; 1,0; 0,0]; elem2 = [1,2,3,4];
for i = 1:refine, [node1, elem1] = uniformrefinequad(node1, elem1); [node2, elem2] = uniformrefinequad(node2, elem2); end

area1 = 1/size(elem1,1) * ones(size(elem1,1),1); area2 = 1/size(elem2,1) * ones(size(elem2,1),1);
[node1, elem1] = generatePerturbingQuadrilateralMesh(node1, elem1, area1, 0.4); 
[node2, elem2] = generatePerturbingQuadrilateralMesh(node2, elem2, area2, 0.4); 

% figure; h = patch('Faces', elem1, 'Vertices', node1); set(h,'facecolor','w','edgecolor','k');
% hold on; h = patch('Faces', elem2, 'Vertices', node2); set(h,'facecolor','w','edgecolor','k');
% view(2); axis equal; axis tight; axis off; hold off;

node = [node1; node2]; elem = [elem1; elem2+size(node1,1)];
dim = ones(1,size(elem,1)); elem = mat2cell(elem, dim); 
weakGalerkin_P1P0P0_polygon(node, elem);
clear all

%% Non-convex mesh
% refine = 3;
% [node1, elem1] = squarequadmesh([0,1,0,1], 1/2^refine); [node1, elem1] = non_convex_octagona_mesh(node1,elem1);
% [node2, elem2] = squarequadmesh([0,1,-1,0], 1/2^refine); [node2, elem2] = non_convex_octagona_mesh(node2,elem2);
% 
% figure; PolyMshr_PlotMsh(node1, elem1, length(elem1)); 
% hold on; PolyMshr_PlotMsh(node2, elem2, length(elem2)); view(2); axis equal; axis tight; axis off; hold off;
% 
% 
% %% Figure of mesh
%     function PolyMshr_PlotMsh(Node,Element,NElem,Supp,Load)
%         % clf; axis equal; axis off; hold on;
%         % Element = Element(1:NElem)';                 %Only plot the first block
%         MaxNVer = max(cellfun(@numel,Element));      %Max. num. of vertices in mesh
%         PadWNaN = @(E) [E' NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
%         ElemMat = cellfun(PadWNaN,Element,'UniformOutput',false);
%         ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
%         
%         patch('Faces',ElemMat,'Vertices',Node,'FaceColor','w'); pause(1e-6)
%         if exist('Supp','var')&&~isempty(Supp) %Plot Supp BC if specified
%             plot(Node(Supp(:,1),1),Node(Supp(:,1),2),'b>','MarkerSize',8);
%         end
%         if exist('Load','var')&&~isempty(Load) %Plot Load BC if specified
%             plot(Node(Load(:,1),1),Node(Load(:,1),2),'m^','MarkerSize',8);
%         end
%     end
end