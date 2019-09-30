function [f1h, f2h] = modifyLoadVectorForLaplaceTerm(f1h, f2h, uLef, uRig, vTop, vBot, uTop, uBot, vLef, vRig, nu)

n = size(f1h,1); h = 1/n;

f1h(:,[1 end]) = [uLef  uRig]; f2h([1 end],:) = [vTop;  vBot]; 

f1h(:,[2 end-1]) = f1h(:,[2 end-1]) + nu/h^2*[uLef  uRig]; 
f2h([2 end-1],:) = f2h([2 end-1],:) + nu/h^2*[vTop; vBot]; 

j = 2:n; f1h(1,j) = f1h(1,j) + nu/h^2*2*uTop(1,j); f1h(end,j) = f1h(end,j) + nu/h^2*2*uBot(1,j);
i = 2:n; f2h(i,1) = f2h(i,1) + nu/h^2*2*vLef(i,1); f2h(i,end) = f2h(i,end) + nu/h^2*2*vRig(i,1);

end