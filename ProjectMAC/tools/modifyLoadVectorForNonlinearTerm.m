function [f1h, f2h] = modifyLoadVectorForNonlinearTerm(f1h, f2h, uLef, uRig, uTop, uBot, vLef, vRig, vTop, vBot, au, bu, av, bv)

n = size(au,1); h = 1/n;

%% au
i = 1:n; 
f1h(i,2) = f1h(i,2) + au(i,2)/(2*h).*uLef(i,1);
f1h(i,n) = f1h(i,n) - au(i,n)/(2*h).*uRig(i,1);

%% bu
j = 2:n;
f1h(1,j) = f1h(1,j) - 2*bu(1,j)/(2*h).*uTop(1,j);
f1h(n,j) = f1h(n,j) + 2*bu(n,j)/(2*h).*uBot(1,j);

%% av
i = 2:n;
f2h(i,1) = f2h(i,1) + 2*av(i,1)/(2*h).*vLef(i,1);
f2h(i,n) = f2h(i,n) - 2*av(i,n)/(2*h).*vRig(i,1);

%% bv
j = 1:n;
f2h(2,j) = f2h(2,j) - bv(2,j)/(2*h).*vTop(1,j);
f2h(n,j) = f2h(n,j) + bv(n,j)/(2*h).*vBot(1,j);


end