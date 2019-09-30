function A = formSp(node, elem)

NT = size(elem,1);
ii = zeros(9*NT,1);  jj = zeros(9*NT,1);  sA = zeros(9*NT,1);
ve(:,:,3) = node(elem(:,2),:)-node(elem(:,1),:);
ve(:,:,1) = node(elem(:,3),:)-node(elem(:,2),:);
ve(:,:,2) = node(elem(:,1),:)-node(elem(:,3),:);
area = 0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));
 
index = 0;
for i = 1:3
    for j = 1:3
        ii(index+1:index+NT) = elem(:,j);
        jj(index+1:index+NT) = elem(:,i);
        sA(index+1:index+NT) = dot(ve(:,:,j),ve(:,:,i),2)./(4*area);
        index = index+NT;
    end
end
N = size(node,1);
A = sparse(ii,jj,sA,N,N);

clear NT ii jj sA ve area 

end