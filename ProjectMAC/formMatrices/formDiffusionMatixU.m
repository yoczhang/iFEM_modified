function DA = formDiffusionMatixU(nu, N, width)
% FORMDIFFUSIONMATRIX Form diffsuion Matrix for NS eqns.
 
%% Parameters and classify fixed dof.
Ndof = N*(N+1); 
h = width/N;
r = nu/h^2;
uidxmat = reshape(1:Ndof,N,N+1);
fixedDof = uidxmat(:,[1 end]);
nearbdDof = uidxmat([1 end],2:end-1);

%% Construct the diffusion matrix
% Matrix for freedof which is a sparse matrix with width 5.
a0 = 4*r*ones(N,N+1);
a0(fixedDof(:)) = 1;     % boundary nodes: the block matrix is I
% a0(nearbdDof(:)) = 5*r;
a0(nearbdDof(:)) = 5*r;

l1 = -r*ones(N,N+1);
l1(fixedDof(:)) = 0;
l1(N,:) = 0;  % l1 is the coefficient of (k,k-1). N+1 and N is not connected

r1 = -r*ones(N,N+1);
r1(fixedDof(:)) = 0;
r1(1,:) = 0; % r1 is the coefficient of (k,k+1). N and N+1 is not connected

l2 = -r*ones(N,N+1);
l2(:,[1 end-1 end]) = 0; % stencil to the left

r2 = -r*ones(N,N+1);
r2(:,[1 2 end]) = 0;     % stencil to the right

DA = spdiags([l2(:),l1(:),a0(:),r1(:),r2(:)],[-N -1:1 N],Ndof,Ndof);