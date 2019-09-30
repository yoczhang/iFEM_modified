function DA = formDiffusionMatixV(nu, N, width)
% FORMDIFFUSIONMATRIX Form diffsuion Matrix for NS eqns. DA = 
 
%% Parameters and classify fixed dof.
Ndof = N*(N+1); 
h = width/N;
r = nu/h^2;
uidxmat = reshape(1:Ndof,N+1,N);
fixedDof = uidxmat([1 end],:);
nearbdDof = uidxmat(2:end-1,[1 end]);

%% Construct the diffusion matrix
% Matrix for freedof which is a sparse matrix with width 5.
a0 = 4*r*ones(N+1,N);
a0(fixedDof(:)) = 1;     % boundary nodes: the block matrix is I
% a0(nearbdDof(:)) = 5*r;
a0(nearbdDof(:)) = 5*r;
l1 = -r*ones(N+1,N);
l1(fixedDof(:)) = 0;
l1(:,N) = 0;  % l1 is the coefficient of (k,k-1). N+1 and N is not connected
r1 = -r*ones(N+1,N);
r1(fixedDof(:)) = 0;
r1(:,1) = 0; % r1 is the coefficient of (k,k+1). N and N+1 is not connected
l2 = -r*ones(N+1,N);
l2([1 end-1 end],:) = 0; % stencil to the left
r2 = -r*ones(N+1,N);
r2([1 2 end],:) = 0;     % stencil to the right
DA = spdiags([l2(:),l1(:),a0(:),r1(:),r2(:)],[-1 -(N+1) 0 N+1 1],Ndof,Ndof);

