function [A, B, C, area] = assembleStiffnessMatrixForPrimalDarcyP0P0RT0(node, elem, elem2edge, NT, NE)

% [A  B;
%  B' C]
% A = sparse(NT,NT); B = sparse(NT,NE); 
C = sparse(NE,NE);

%% 
% for n = 1:NT
%     vertices = node(elem(n,:),:); % 3*2
%     area = abs(det([1 1 1; vertices']))/2;
%     e1 = sqrt(sum((node(elem(n,2),:) - node(elem(n,3),:)).^2, 2));   
%     e2 = sqrt(sum((node(elem(n,3),:) - node(elem(n,1),:)).^2, 2));   
%     e3 = sqrt(sum((node(elem(n,1),:) - node(elem(n,2),:)).^2, 2));  
%     l1 = e1^2; l2 = e2^2; l3 = e3^2;
%     l12 = e1^2 + e2^2; l23 = e2^2 + e3^2; l13 = e3^2 + e1^2;
%     l123 = e1^2 + e2^2 + e3^2;
%     M00 = 144*area/l123;
%     Mob = -48*area/l123*ones(1,3);
%     Mbb = 16*area/l123*ones(3,3) + ...
%           1/(2*area)*[2*l1 l3-l12 l2-l13; l3-l12 2*l2 l1-l23; l2-l13 l1-l23 2*l3];
%     
%     A(n,n) = A(n,n) + M00;
%     for alpha = 1:3
%         B(n,elem2edge(n,alpha)) = B(n,elem2edge(n,alpha)) + Mob(1,alpha);
%         for beta = 1:3
%             C(elem2edge(n,beta),elem2edge(n,alpha)) = ...
%                 C(elem2edge(n,beta),elem2edge(n,alpha)) + ...
%                 Mbb(alpha,beta);
%         end
%     end
% end


%%
% % compute ct2 = 1/mean(||x-xc||^2)
% center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
% mid1 = (node(elem(:,2),:) + node(elem(:,3),:))/2;
% mid2 = (node(elem(:,3),:) + node(elem(:,1),:))/2;
% mid3 = (node(elem(:,1),:) + node(elem(:,2),:))/2;
% ct2 = 3./sum((mid1 - center).^2 + (mid2 - center).^2 + (mid3 - center).^2,2);
% [Dphi,area] = gradbasis(node,elem);
% clear center mid1 mid2 mid3
% 
% % Mbb: edge - edge           
% for i = 1:3
%     for j = i:3
%         % local to global index map
%         ii = double(elem2edge(:,i));
%         jj = double(elem2edge(:,j));
%         % local stiffness matrix
%         Cij = 4*dot(Dphi(:,:,i),Dphi(:,:,j),2).*area + 4/9*ct2.*area;      
%         if (j==i)
%             C = C + sparse(ii, jj, Cij, NE, NE);
%         else
%             C = C + sparse([ii,jj], [jj,ii], [Cij; Cij], NE,NE);        
%         end        
%     end
% end
% 
% % Mob: interior - edge
% Bij = -4/3*ct2.*area;
% B = sparse([(1:NT)', (1:NT)', (1:NT)'], ...
%              double(elem2edge(:)), [Bij, Bij, Bij], NT, NE);
% 
% % Moo: diagonal of interor
% Aij = 4*ct2.*area;
% A = sparse(1:NT, 1:NT, Aij, NT, NT);
% 
% clear Dphi area Cij Bij Aij

%%
ve1 = node(elem(:,3),:)-node(elem(:,2),:);
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
area = 0.5*(-ve3(:,1).*ve2(:,2) + ve3(:,2).*ve2(:,1));
area = abs(area);

l1 = sum((node(elem(:,3),:) - node(elem(:,2),:)).^2, 2);
l2 = sum((node(elem(:,1),:) - node(elem(:,3),:)).^2, 2);
l3 = sum((node(elem(:,2),:) - node(elem(:,1),:)).^2, 2);
l12 = l1 + l2;
l23 = l2 + l3;
l31 = l3 + l1; l13 = l31;
l123 = l1 + l2 + l3;

A = spdiags(144*area./l123, 0, NT, NT);

B = sparse([(1:NT)', (1:NT)', (1:NT)'], ...
             double(elem2edge(:)), [-48*area./l123, -48*area./l123, -48*area./l123], NT, NE);

Mbb = cell(3,3);
Mbb{1,1} = 16*area./l123 + 1./(2*area)*2.*l1;
Mbb{1,2} = 16*area./l123 + 1./(2*area).*(l3 - l12);
Mbb{1,3} = 16*area./l123 + 1./(2*area).*(l2 - l13);
Mbb{2,2} = 16*area./l123 + 1./(2*area)*2.*l2;
Mbb{2,3} = 16*area./l123 + 1./(2*area).*(l1 - l23);
Mbb{3,3} = 16*area./l123 + 1./(2*area)*2.*l3;
Mbb{2,1} = Mbb{1,2}; Mbb{3,1} = Mbb{1,3}; Mbb{3,2} = Mbb{2,3};
for i = 1:3
    for j = 1:3
        % local to global index map
        ii = double(elem2edge(:,i));
        jj = double(elem2edge(:,j));
        % local stiffness matrix
        Cij = Mbb{i,j};      
%         if (j==i)
            C = C + sparse(ii, jj, Cij, NE, NE);
%         else
%             C = C + sparse([ii,jj], [jj,ii], [Cij; Cij], NE,NE);        
%         end        
    end
end
end