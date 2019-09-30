function [au, bu, av, bv] = formConvectionCoefficientMatrices(uh, vh, uexact, vexact, level)

%    [ au   bu;
%      av   bv ] 
n = size(uh,1);

%% For direct solver
% h = 1 / n; 
% [ux,uy] = meshgrid(0: h: 1, 1-h/2: -h: h/2);
% [vx,vy] = meshgrid(h/2: h: 1-h/2, 1: -h: 0);
% au = uexact(ux(:), uy(:)); bu = vexact(ux(:), uy(:));
% av = uexact(vx(:), vy(:)); bv = vexact(vx(:), vy(:));
% au = reshape(au,n,n+1); bu = reshape(bu,n,n+1);
% av = reshape(av,n+1,n); bv = reshape(bv,n+1,n);
% 
% clear ux uy vx vy
%%
au = cell(level,1); av = cell(level,1); bu = cell(level,1); bv = cell(level,1);
for j = 1:level
    n = 2^(j+1); h = 1 / n; 
    [ux,uy] = meshgrid(0: h: 1, 1-h/2: -h: h/2);
    [vx,vy] = meshgrid(h/2: h: 1-h/2, 1: -h: 0);
    au{j} = uexact(ux(:), uy(:)); bu{j} = vexact(ux(:), uy(:));
    av{j} = uexact(vx(:), vy(:)); bv{j} = vexact(vx(:), vy(:));
    au{j} = reshape(au{j},n,n+1); bu{j} = reshape(bu{j},n,n+1);
    av{j} = reshape(av{j},n+1,n); bv{j} = reshape(bv{j},n+1,n);
end
clear ux uy vx vy

%% For wcycle
% au = cell(level,1); av = cell(level,1); bu = cell(level,1); bv = cell(level,1);
% 
% %% The level_th
% au{level} = uh;  bv{level} = vh;
% 
% n = 2^(level+1); h = 1/n; bu{level} = zeros(n,n+1); av{level} = zeros(n+1,n);
% 
% i = 1:n; j = 2:n; bu{level}(i,j) = (bv{level}(i,j-1) + bv{level}(i+1,j-1) + bv{level}(i,j) + bv{level}(i+1,j)) / 4;
% i = 2:n; j = 1:n; av{level}(i,j) = (au{level}(i-1,j) + au{level}(i,j) + au{level}(i-1,j+1) + au{level}(i,j+1)) / 4;
% 
% yy = (1-h/2:-h:h/2)'; vLef = vexact(0,yy); vRig = vexact(1,yy); bu{level}(:,[1,end]) = [vLef vRig];
% xx = h/2:h:1-h/2;  i = 1:n; uTop(1,i) = uexact(xx(i),1); uBot(1,i) = uexact(xx(i),0); av{level}([1,end],:) = [uTop; uBot];
% 
% clear n h i j xx yy vLef vRig uTop uBot
% %% from (level-1)_th to 1th
% for k = level-1:-1:1
%     n = 2^(k+1); h = 1/n;
%     [au{k}, bv{k}] = Res(au{k+1}, bv{k+1}, sparse(n^2,n^2));
%     
%     % impose Dirichlet boundary condition
%     yy = (1-h/2:-h:h/2)'; uLef = uexact(0,yy); uRig = uexact(1,yy); au{k}(:,[1,end]) = [uLef uRig];
%     xx = h/2:h:1-h/2;  vTop = vexact(xx,1); vBot = vexact(xx,0); bv{k}([1,end],:) = [vTop; vBot];
%     
%     % form bu and av
%     i = 1:n; j = 2:n; bu{k}(i,j) = (bv{k}(i,j-1) + bv{k}(i+1,j-1) + bv{k}(i,j) + bv{k}(i+1,j)) / 4;
%     i = 2:n; j = 1:n; av{k}(i,j) = (au{k}(i-1,j) + au{k}(i,j) + au{k}(i-1,j+1) + au{k}(i,j+1)) / 4;
%     yy = (1-h/2:-h:h/2)'; vLef = vexact(0,yy); vRig = vexact(1,yy); bu{k}(:,[1,end]) = [vLef vRig];
%     xx = h/2:h:1-h/2;  i = 1:n; uTop(1,i) = uexact(xx(i),1); uBot(1,i) = uexact(xx(i),0); av{k}([1,end],:) = [uTop; uBot];
%     clear n h yy xx uLef uRig vTop vBot vLef vRig uTop uBot 
% end

end