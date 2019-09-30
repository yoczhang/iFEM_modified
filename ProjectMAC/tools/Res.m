function [rH1, rH2, rH3] = Res(rh1, rh2, rh3)

N  = size(rh1,1) / 2;
rH1 = zeros(N,N+1); rH2 = zeros(N+1,N); rH3 = zeros(N,N);

i = 1:N; j = 2:N;
rH1(i,j) = (rh1(2*i-1,2*j-2) + rh1(2*i,2*j-2) + 2*rh1(2*i-1,2*j-1) + 2*rh1(2*i,2*j-1) + rh1(2*i-1,2*j) + rh1(2*i,2*j))/8;

i = 2:N; j = 1:N;
rH2(i,j) = (rh2(2*i-2,2*j-1) + 2*rh2(2*i-1,2*j-1) + rh2(2*i,2*j-1) + rh2(2*i-2,2*j) + 2*rh2(2*i-1,2*j) + rh2(2*i,2*j))/8;

i = 1:N; j = 1:N;
rH3(i,j) = (rh3(2*i-1,2*j-1) + rh3(2*i,2*j-1) + rh3(2*i-1,2*j) + rh3(2*i,2*j))/4;  
        
end