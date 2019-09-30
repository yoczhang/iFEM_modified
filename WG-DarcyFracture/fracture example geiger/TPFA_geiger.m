function [r, b] = TPFA_geiger(NGamma, h, n, lGamma, KGammaTau)

co = lGamma*KGammaTau; alpha = 2*co/h;
r = sparse(NGamma, NGamma); b = zeros(NGamma,1);

N = 2^n;
boun = [1, 4*N, 4*N+1, 5*N, 5*N+1, 6*N, 6*N+1, 8*N, ...
        8*N+1, 12*N, 12*N+1, 13*N, 13*N+1, 14*N, 14*N+1, 16*N, ...
        16*N+1, 17*N, 17*N+1, 18*N, 18*N+1, 20*N, ...
        20*N+1, 21*N, 21*N+1, 22*N, 22*N+1, 24*N, ...
        24*N+1, 25*N, 25*N+1, 26*N, ...
        26*N+1, 27*N, 27*N+1, 28*N];
inte = setdiff(1:NGamma, boun);

%% 统一处理内部单元
for i = 1:length(inte)
    j = inte(i);
    r(j,j) = r(j,j) + alpha; r(j,j-1) = r(j,j-1) - alpha/2; r(j,j+1) = r(j,j+1) - alpha/2;
end

%% 左侧Neuamnn边界, 边界值为 1
r(1,1) = r(1,1) + alpha/2; r(1,2) = r(1,2) - alpha/2; 
% b(1) = 1;

%% 上下Neuamnn边界， 边界值为 0
r(16*N,16*N)   = r(16*N,16*N)   + alpha/2; r(16*N,16*N-1) = r(16*N,16*N-1) - alpha/2;
r(8*N+1,8*N+1) = r(8*N+1,8*N+1) + alpha/2; r(8*N+1,8*N+2) = r(8*N+1,8*N+2) - alpha/2;
r(24*N,24*N)   = r(24*N,24*N)   + alpha/2; r(24*N,24*N-1) = r(24*N,24*N-1) - alpha/2;

%% 右侧Dirichlet边界， 边界值为1
r(8*N,8*N)   = r(8*N,8*N)   + 3*alpha/2; r(8*N,8*N-1)   = r(8*N,8*N-1)   - alpha/2; b(8*N)  = alpha;
r(20*N,20*N) = r(20*N,20*N) + 3*alpha/2; r(20*N,20*N-1) = r(20*N,20*N-1) - alpha/2; b(20*N) = alpha;

%% 处理内部交点
angle1 = [0     pi/2  0    pi/2;
          pi/2  0     pi/2  0;
          0    pi/2  0     pi/2;
          pi/2  0    pi/2  0];
% angle1 = [0     pi/2  pi    pi/2;
%           pi/2  0     pi/2  pi;
%           pi    pi/2  0     pi/2;
%           pi/2  pi    pi/2  0];
inter1 = [2  10 3  11;
          30 34 31 35;
          20 26 21 27];
for i = 1:size(inter1,1)
    j = boun(inter1(i,:));
    for k = 1:4
        temp = setdiff(1:4, k);
        r(j(k),j(k)) = r(j(k),j(k)) + alpha/4*(cos(angle1(k,temp(1))/2) + cos(angle1(k,temp(2))/2) + cos(angle1(k,temp(3))/2));
        r(j(k),j(temp(1))) = r(j(k),j(temp(1))) - alpha/4*cos(angle1(k,temp(1))/2);
        r(j(k),j(temp(2))) = r(j(k),j(temp(2))) - alpha/4*cos(angle1(k,temp(2))/2);
        r(j(k),j(temp(3))) = r(j(k),j(temp(3))) - alpha/4*cos(angle1(k,temp(3))/2);
    end
end


angle2 = [0     pi/2  0;
          pi/2  0     pi/2;
          0    pi/2  0];
% angle2 = [0     pi/2  pi;
%           pi/2  0     pi/2;
%           pi    pi/2  0];
inter2 = [13 29 12;
          15 17 14;
          5  33 4;
          19 36 18;
          7  23 6;
          25 32 24];
for i = 1:size(inter2,1)
    j = boun(inter2(i,:));
    for k = 1:3
        temp = setdiff(1:3, k);
        r(j(k),j(k)) = r(j(k),j(k)) + alpha/3*(cos(angle2(k,temp(1))/2) + cos(angle2(k,temp(2))/2));
        r(j(k),j(temp(1))) = r(j(k),j(temp(1))) - alpha/3*cos(angle2(k,temp(1))/2);
        r(j(k),j(temp(2))) = r(j(k),j(temp(2))) - alpha/3*cos(angle2(k,temp(2))/2);
    end
end

i = [2, 4, 6, 10, 12, 14, 18, 20, 24, 26, 30, 32, 34, 36];
j = [3, 5, 7, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35];
i = boun(i)'; j = boun(j)';
ii = [i; i;   j; j];
jj = [i; i-1; j; j+1];
ss = [  alpha/2*ones(length(i),1);
      - alpha/2*ones(length(i),1);
        alpha/2*ones(length(j),1);
      - alpha/2*ones(length(j),1)];
r = r + sparse(ii, jj, ss, NGamma, NGamma);

end