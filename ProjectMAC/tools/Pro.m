function [duh, dvh, dph] = Pro(duH, dvH, dpH)

N = size(duH,1); n = 2*N;

duh = zeros(n,n+1); dvh = zeros(n+1,n); dph = zeros(n,n);

%% u
i = 1:N-1; j = 2:N;
duh(2*i,2*j-1) = 0.75*duH(i,j) + 0.25*duH(i+1,j);
duh(2*i+1,2*j-1) = 0.25*duH(i,j) + 0.75*duH(i+1,j);

i = 1:N-1; j = 2:N+1;
duh(2*i,2*j-2) = 0.5*duh(2*i,2*j-3) + 0.5*duh(2*i,2*j-1);
duh(2*i+1,2*j-2) = 0.5*duh(2*i+1,2*j-3) + 0.5*duh(2*i+1,2*j-1);

j = 2:N; duh(1,2*j-1) = 0.5*duH(1,j); duh(n,2*j-1) = 0.5*duH(end,j); 
j = 2:N+1; duh(1,2*j-2) = 0.5*duh(1,2*j-3) + 0.5*duh(1,2*j-1);
duh(n,2*j-2) = 0.5*duh(n,2*j-3) + 0.5*duh(n,2*j-1);

% i = 1:N-1; duh(2*i,1) = 0.75*duH(i,1) + 0.25*duH(i+1,1); duh(2*i+1,1) = 0.25*duH(i,1) + 0.75*duH(i+1,1);
% duh(1,1) = 0.5*duH(1,1); duh(n,1) = 0.5*duH(N,1);
% i = 1:N-1; duh(2*i,n+1) = 0.75*duH(i,N+1) + 0.25*duH(i+1,N+1); duh(2*i+1,n+1) = 0.25*duH(i,N+1) + 0.75*duH(i+1,N+1);
% duh(1,n+1) = 0.5*duH(1,N+1); duh(n,n+1) = 0.5*duH(N,N+1);


%% v
i = 2:N; j = 1:N-1;
dvh(2*i-1,2*j) = 0.75*dvH(i,j) + 0.25*dvH(i,j+1);
dvh(2*i-1,2*j+1) = 0.25*dvH(i,j) + 0.75*dvH(i,j+1);

i = 2:N+1; j = 1:N-1;
dvh(2*i-2,2*j) = 0.5*dvh(2*i-3,2*j) + 0.5*dvh(2*i-1,2*j);
dvh(2*i-2,2*j+1) = 0.5*dvh(2*i-3,2*j+1) + 0.5*dvh(2*i-1,2*j+1);

i = 2:N; dvh(2*i-1,1) = 0.5*dvH(i,1); dvh(2*i-1,n) = 0.5*dvH(i,end); 
i = 2:N+1; dvh(2*i-2,1) = 0.5*dvh(2*i-3,1) + 0.5*dvh(2*i-1,1);
dvh(2*i-2,n) = 0.5*dvh(2*i-3,n) + 0.5*dvh(2*i-1,n);


% i = 1:N-1; dvh(1,2*i) = 0.75*dvH(1,i) + 0.25*dvH(1,i+1); dvh(1,2*i+1) = 0.25*dvH(1,i) + 0.75*dvH(1,i+1);
% dvh(1,1) = 0.5*dvH(1,1); dvh(1,n) = 0.5*dvH(1,N);
% i = 1:N-1; dvh(n+1,2*i) = 0.75*dvH(N+1,i) + 0.25*dvH(N+1,i+1); dvh(n+1,2*i+1) = 0.25*dvH(N+1,i) + 0.75*dvH(N+1,i+1);
% dvh(n+1,1) = 0.5*dvH(N+1,1); dvh(n+1,n) = 0.5*dvH(N+1,N);


%% p
i = 1:N; j = 1:N;
dph(2*i-1,2*j-1) = dpH(i,j); dph(2*i,2*j-1) = dpH(i,j);
dph(2*i-1,2*j) = dpH(i,j); dph(2*i,2*j) = dpH(i,j);
end