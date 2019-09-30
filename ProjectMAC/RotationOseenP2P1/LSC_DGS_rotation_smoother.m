function r = LSC_DGS_rotation_smoother(bigu, A, B, bigF, Sp, node, elem)

Nu = size(A,1)/2; Np = size(B,1); i = 1:Np-1;
A1 = A(1:Nu,1:Nu); M = A(Nu+1:2*Nu,1:Nu);
u = bigu(1:2*Nu,1); p = bigu(2*Nu+1:end,1);

ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
area = 0.5*(-ve3(:,1).*ve2(:,2) + ve3(:,2).*ve2(:,1));
clear ve2 ve3

%% Step 1, Relax the momentum equation
B1 = B(:,1:Nu); B2 = B(:,Nu+1:2*Nu);
F = bigF(1:2*Nu,1) - A*u - [B1'*p; B2'*p];
u = u + A \ F;


%% Step 2, Relax the transformed continuity equation
g = bigF(2*Nu+1:end,1);
rq = g - B*u;
dq = zeros(Np,1);
% dq = dq + tril(Sp) \ (rq - Sp*dq);
% dq = dq + triu(Sp) \ (rq - Sp*dq);
dq(i) = Sp(i,i) \ rq(i);
c = sum(mean(dq(elem),2).*area)/sum(area);
dq = dq - c;

%% Step 3, Distribute the correction back to the original variables
u = u + B'*dq;

rq2 = B*A*B'*dq;
ep = zeros(Np,1);
% ep = ep + tril(Sp) \ (rq2 - Sp*ep);
% ep = ep + triu(Sp) \ (rq2 - Sp*ep);
ep(i) = Sp(i,i) \ rq2(i);
c = sum(mean(ep(elem),2).*area)/sum(area);
ep = ep - c;
p = p - ep;

r = [u; p];

end