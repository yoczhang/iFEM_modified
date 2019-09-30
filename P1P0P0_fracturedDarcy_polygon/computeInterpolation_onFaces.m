function r = computeInterpolation_onFaces(node,edge,NE,p1,p2)

t = [0, sqrt(21)/7, -sqrt(21)/7, 1, -1];  weight = [32/45, 49/90, 49/90, 1/10, 1/10];

%% 
r1 = zeros(NE/2,1); r2 = zeros(NE/2,1); left = node(edge(:,1),:); right = node(edge(:,2),:);

for k = 1:length(weight)
    pxy1 = (right(1:NE/2,:) - left(1:NE/2,:))/2 * t(k) + (right(1:NE/2,:) + left(1:NE/2,:))/2;
    pxy2 = (right(NE/2+1:end,:) - left(NE/2+1:end,:))/2 * t(k) + (right(NE/2+1:end,:) + left(NE/2+1:end,:))/2;
    
    r1 = r1 + 1/2 .* weight(k) .* p1(pxy1); r2 = r2 + 1/2 .* weight(k) .* p2(pxy2);
end

r = [r1; r2];

%%
% r = zeros(NE,1); 
% left = node(edge(:,1),:); right = node(edge(:,2),:);
% 
% for k = 1:length(weight)
%     pxy1 = (right - left)/2 * t(k) + (right + left)/2;
%     
%     r = r + 1/2 .* weight(k) .* p1(pxy1);
% end

end