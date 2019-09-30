function [u, v] = symmetricGaussSeidel_Velocity(solu, A, B, t, dofu)

n = sqrt(size(B,1));

uv = solu(1:2*dofu); p = solu(2*dofu+1:end);
bigF = t(1:2*dofu) - B'*p;

uv = uv + tril(A) \ (bigF - A*uv);
uv = uv + triu(A) \ (bigF - A*uv);

u = reshape(uv(1:dofu),n,n+1);
v = reshape(uv(dofu+1:2*dofu),n+1,n);

%% Forward for u
% for j = 2:n
%     u(1,j) = (infFlow*(2*uTop(1,j) + u(2,j) + u(1,j-1) + u(1,j+1)) - au(1,j).*(u(1,j+1) - u(1,j-1)) - ...
%               bu(1,j).*(2*uTop(1,j) - u(2,j)) - 2*(p(1,j) - p(1,j-1)) + 2*h*f1h(1,j) ) / (5*infFlow - bu(1,j));
%     for i = 2:n-1
%         u(i,j) = (2*(h^2*f1h(i,j) - h*(p(i,j) - p(i,j-1)) - 0.5*h*au(i,j).*(u(i,j+1)- u(i,j-1)) - ...
%                  0.5*h*bu(i,j).*(u(i-1,j) - u(i+1,j)) ) / (h*infFlow) + ...
%                   u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1)) / 4;
%     end
%     u(n,j) = (infFlow*(u(n-1,j) + 2*uBot(1,j) + u(n,j-1) + u(n,j+1)) - au(n,j).*(u(n,j+1) - u(n,j-1)) - ...
%               bu(n,j).*(u(n-1,j) - 2*uBot(1,j)) - 2*(p(n,j) - p(n,j-1)) + 2*h*f1h(n,j) ) / (5*infFlow + bu(n,j));
% end
% 
% 
% %% Backward for u
% for j = n:-1:2
%     u(n,j) = (infFlow*(u(n-1,j) + 2*uBot(1,j) + u(n,j-1) + u(n,j+1)) - au(n,j).*(u(n,j+1) - u(n,j-1)) - ...
%               bu(n,j).*(u(n-1,j) - 2*uBot(1,j)) - 2*(p(n,j) - p(n,j-1)) + 2*h*f1h(n,j) ) / (5*infFlow + bu(n,j));  
%     for i = n-1:-1:2
%         u(i,j) = (2*(h^2*f1h(i,j) - h*(p(i,j) - p(i,j-1)) - 0.5*h*au(i,j).*(u(i,j+1)- u(i,j-1)) - ...
%                   0.5*h*bu(i,j).*(u(i-1,j) - u(i+1,j)) ) / (h*infFlow) + ...
%                   u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1)) / 4;
%     end
%     u(1,j) = (infFlow*(2*uTop(1,j) + u(2,j) + u(1,j-1) + u(1,j+1)) - au(1,j).*(u(1,j+1) - u(1,j-1)) - ...
%               bu(1,j).*(2*uTop(1,j) - u(2,j)) - 2*(p(1,j) - p(1,j-1)) + 2*h*f1h(1,j) ) / (5*infFlow - bu(1,j)); 
% end
% 
% 
%% Forward for v
% for i = 2:n
%         v(i,1) = (infFlow*(v(i-1,1) + v(i+1,1) + 2*vLef(i,1) + v(i,2)) - av(i,1).*(v(i,2) - 2*vLef(i,1)) - ...
%                   bv(i,1).*(v(i-1,1) - v(i+1,1)) - 2*(p(i-1,1) - p(i,1)) + 2*h*f2h(i,1)) / (5*infFlow + av(i,1));
% end
% for j = 2:n-1  
%     for i = 2:n
%         v(i,j) = (2*(h^2*f2h(i,j) - h*(p(i-1,j) - p(i,j)) - 0.5*h*bv(i,j).*(v(i-1,j)- v(i+1,j)) - ...
%                   0.5*h*av(i,j).*(v(i,j+1)- v(i,j-1))) / (h*infFlow) + ...
%                   v(i-1,j) + v(i+1,j) + v(i,j-1) + v(i,j+1)) / 4;
%     end
% end
% for i = 2:n
%     v(i,n) = (infFlow*(v(i-1,n) + v(i+1,n) + v(i,n-1) + 2*vRig(i,1)) - av(i,n).*(2*vRig(i,1) - v(i,n-1)) - ...
%               bv(i,n).*(v(i-1,n) - v(i+1,n)) - 2*(p(i-1,n) - p(i,n)) + 2*h*f2h(i,n)) / (5*infFlow - av(i,n));
% end


%% Backward for v
% for i = n:-1:2
%     v(i,n) = (infFlow*(v(i-1,n) + v(i+1,n) + v(i,n-1) + 2*vRig(i,1)) - av(i,n).*(2*vRig(i,1) - v(i,n-1)) - ...
%               bv(i,n).*(v(i-1,n) - v(i+1,n)) - 2*(p(i-1,n) - p(i,n)) + 2*h*f2h(i,n)) / (5*infFlow - av(i,n));      
% end
% for j = n-1:-1:2  
%     for i = n:-1:2
%         v(i,j) = (2*(h^2*f2h(i,j) - h*(p(i-1,j) - p(i,j)) - 0.5*h*bv(i,j).*(v(i-1,j)- v(i+1,j)) - ...
%                   0.5*h*av(i,j).*(v(i,j+1)- v(i,j-1))) / (h*infFlow) + ...
%                   v(i-1,j) + v(i+1,j) + v(i,j-1) + v(i,j+1)) / 4;
%     end
% end
% for i = n:-1:2
%         v(i,1) = (infFlow*(v(i-1,1) + v(i+1,1) + 2*vLef(i,1) + v(i,2)) - av(i,1).*(v(i,2) - 2*vLef(i,1)) - ...
%                   bv(i,1).*(v(i-1,1) - v(i+1,1)) - 2*(p(i-1,1) - p(i,1)) + 2*h*f2h(i,1)) / (5*infFlow + av(i,1));
% end


end