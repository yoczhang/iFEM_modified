function [a1, a2] = formConvectionMatrices(au, bu, av, bv)

n = size(au,1); dof = n*(n+1); h = 1/n; 

a = 1:dof; a = reshape(a,n,n+1);
%% au
ii1 = a(:,2); jj1 = a(:,3); ss1 = au(:,2)/(2*h); % lef
ii2 = a(:,n); jj2 = a(:,n-1); ss2 = - au(:,n)/(2*h); % rig

i = 1:n; j = 3:n-1; ii3 = a(i,j); jj3Lef = a(i,j-1); jj3Rig = a(i,j+1); ss3Lef = - au(i,j)/(2*h); ss3Rig = au(i,j)/(2*h);

r1 = sparse([ii1;ii2;ii3(:);ii3(:)], [jj1;jj2;jj3Lef(:);jj3Rig(:)], [ss1; ss2;ss3Lef(:);ss3Rig(:)], dof, dof);
clear ii1 jj1 ss1 ii2 jj2 ss2 i j ii3 jj3Lef jj3Rig ss3Lef ss3Rig  

%% bu
ii1 = a([1,1],2:n); jj1 = a([1,2],2:n); ss1 = - bu([1,1],2:n)/(2*h); % top
ii2 = a([n,n],2:n); jj2 = a([n-1,n],2:n); ss2 = bu([n,n],2:n)/(2*h); % bot

i = 2:n-1; j = 2:n; ii3 = a(i,j); jj3Top = a(i-1,j); jj3Bot = a(i+1,j); ss3Top = bu(i,j)/(2*h); ss3Bot = - bu(i,j)/(2*h);

r2 = sparse([ii1(:);ii2(:);ii3(:);ii3(:)], [jj1(:);jj2(:);jj3Top(:);jj3Bot(:)], ...
            [ss1(:);ss2(:);ss3Top(:);ss3Bot(:)], dof, dof);
clear ii1 jj1 ss1 ii2 jj2 ss2 i j ii3 jj3Top jj3Bot ss3Top ss3Bot

a = 1:dof; a = reshape(a, n+1, n);
%% av
ii1 = a(2:n,[1,1]); jj1 = a(2:n,[1,2]); ss1 = av(2:n,[1,1])/(2*h); % lef
ii2 = a(2:n,[n,n]); jj2 = a(2:n,[n-1,n]); ss2 = - av(2:n,[n,n])/(2*h); % rig

i = 2:n; j = 2:n-1; ii3 = a(i,j); jj3Lef = a(i,j-1); jj3Rig = a(i,j+1); ss3Lef = - av(i,j)/(2*h); ss3Rig = av(i,j)/(2*h);

r3 = sparse([ii1(:);ii2(:);ii3(:);ii3(:)], [jj1(:);jj2(:);jj3Lef(:);jj3Rig(:)], ...
            [ss1(:);ss2(:);ss3Lef(:);ss3Rig(:)], dof, dof);
clear ii1 jj1 ss1 ii2 jj2 ss2 i j ii3 jj3Lef jj3Rig ss3Lef ss3Rig

%% bv
ii1 = a(2,:); jj1 = a(3,:); ss1 = - bv(2,:)/(2*h); % top
ii2 = a(n,:); jj2 = a(n-1,:); ss2 = bv(n,:)/(2*h); % bot

i = 3:n-1; j = 1:n; ii3 = a(i,j); jj3Top = a(i-1,j); jj3Bot = a(i+1,j); ss3Top = bv(i,j)/(2*h); ss3Bot = - bv(i,j)/(2*h);

r4 = sparse([ii1(:);ii2(:);ii3(:);ii3(:)], [jj1(:);jj2(:);jj3Top(:);jj3Bot(:)], ...
            [ss1(:);ss2(:);ss3Top(:);ss3Bot(:)], dof, dof);
clear ii1 jj1 ss1 ii2 jj2 ss2 i j ii3 jj3Top jj3Bot ss3Top ss3Bot a

%%
a1 = r1 + r2; a2 = r3+ r4;
clear r1 r2 r3 r4
end