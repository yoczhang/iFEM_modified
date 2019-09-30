function [u, v] = symmetricGaussSeidel_Velocity_rotation(u, v, p, f1h, f2h, A1, A2, B, U, V)


n = size(u,1); dofu = n*(n+1); 
u = u(:); v = v(:); vel = [u; v]; p = p(:); 
f1h = f1h(:); f2h = f2h(:); fh = [f1h; f2h];
bigA  = [A1 U; V A2];  B1 = B(:,1:dofu); B2 = B(:,dofu+1:2*dofu);

%% block symmetric Gauss-Seidel smoother
u = u + tril(A1) \ (f1h - A1*u - U*v - B1'*p);
u = u + triu(A1) \ (f1h - A1*u - U*v - B1'*p);

v = v + tril(A2) \ (f2h - A2*v - V*u - B2'*p);
v = v + triu(A2) \ (f2h - A2*v - V*u - B2'*p);

v = v + tril(A2) \ (f2h - A2*v - V*u - B2'*p);
v = v + triu(A2) \ (f2h - A2*v - V*u - B2'*p);

u = u + tril(A1) \ (f1h - A1*u - U*v - B1'*p);
u = u + triu(A1) \ (f1h - A1*u - U*v - B1'*p);

u = reshape(u, n, n+1); v = reshape(v, n+1, n);

% bigF = [f1h - B1'*p; f2h - B2'*p];
% vel = bigA \ bigF;
% u = reshape(vel(1:dofu), n, n+1);
% v = reshape(vel(dofu+1:2*dofu), n+1, n);

%% Collective Symmetric Gauss-Seidel smoother
% smootheru = [triu(A1) triu(U); triu(V) triu(A2)];
% smootherl = [tril(A1) tril(U); tril(V) tril(A2)];
% 
% vel = vel + smootherl \ (fh - bigA*vel - B'*p);
% vel = vel + smootheru \ (fh - bigA*vel - B'*p);
% u = reshape(vel(1:dofu,1), n, n+1);
% v = reshape(vel(dofu+1:2*dofu), n+1, n);


%% Collective Jacobi smoother
% weight = 0.5;
% smoother = [spdiags(diag(A1),0,dofu,dofu)   spdiags(diag(U),0,dofu,dofu);
%             spdiags(diag(V),0,dofu,dofu)    spdiags(diag(A2),0,dofu,dofu)];
% 
% vel = vel + weight*smoother \ (fh - bigA*vel - B'*p);
% u = reshape(vel(1:dofu,1), n, n+1);
% v = reshape(vel(dofu+1:2*dofu), n+1, n);


end