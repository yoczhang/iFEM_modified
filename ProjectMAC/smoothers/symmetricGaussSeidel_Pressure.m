function dq = symmetricGaussSeidel_Pressure(dq, rc, h, Lap)

n = size(rc,1); 

dq = dq(:); rc = rc(:); 
dq = dq + tril(Lap) \ (rc - Lap*dq);
dq = dq + triu(Lap) \ (rc - Lap*dq);

dq = reshape(dq, n, n);

%% Forward and Backward in form of matrix vector

% dq(1,1) = (dq(2,1) + dq(1,2) + h^2*rc(1,1))/2; % left top
% for i = 2:n-1
%     dq(i,1) = (dq(i-1,1) + dq(i+1,1) + dq(i,2) + h^2*rc(i,1))/3; % left 
% end
% dq(n,1) = (dq(n-1,1) + dq(n,2) + h^2*rc(n,1))/2; % left bottom
% for j = 2:n-1
%     dq(1,j) = (dq(2,j) + dq(1,j-1) + dq(1,j+1) + h^2*rc(1,j))/3; % top
%     i = 2:n-1;
%         dq(i,j) = (h^2*rc(i,j) + dq(i-1,j) + dq(i+1,j) + dq(i,j-1) + dq(i,j+1))/4;
% %     end
%     dq(n,j) = (dq(n-1,j) + dq(n,j-1) + dq(n,j+1) + h^2*rc(n,j))/3; % bottom
% end
% dq(1,n) = (dq(2,n) + dq(1,n-1) + h^2*rc(1,n))/2; % right top
% for i = 2:n-1
%     dq(i,n) = (dq(i-1,n) + dq(i+1,n) + dq(i,n-1) + h^2*rc(i,n))/3; % right
% end
% dq(n,n) = (dq(n-1,n) + dq(n,n-1) + h^2*rc(n,n))/2; % right bottom
% 
% 
% % Backward
% 
% dq(n,n) = (dq(n-1,n) + dq(n,n-1) + h^2*rc(n,n))/2; % right bottom
% for i = n-1:-1:2
%     dq(i,n) = (dq(i-1,n) + dq(i+1,n) + dq(i,n-1) + h^2*rc(i,n))/3; % right
% end
% dq(1,n) = (dq(2,n) + dq(1,n-1) + h^2*rc(1,n))/2; % right top
% for j = n-1:-1:2
%     dq(n,j) = (dq(n-1,j) + dq(n,j-1) + dq(n,j+1) + h^2*rc(n,j))/3; % bottom    
%     i = n-1:-1:2;
%         dq(i,j) = (h^2*rc(i,j) + dq(i-1,j) + dq(i+1,j) + dq(i,j-1) + dq(i,j+1))/4;
% %     end
%     dq(1,j) = (dq(2,j) + dq(1,j-1) + dq(1,j+1) + h^2*rc(1,j))/3; % top
% end
% dq(n,1) = (dq(n-1,1) + dq(n,2) + h^2*rc(n,1))/2; % left bottom
% for i = n-1:-1:2
%     dq(i,1) = (dq(i-1,1) + dq(i+1,1) + dq(i,2) + h^2*rc(i,1))/3; % left 
% end
% dq(1,1) = (dq(2,1) + dq(1,2) + h^2*rc(1,1))/2; % left top

end