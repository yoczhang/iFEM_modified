function [rh1, rh2, rh3] = getResidualOseen(uh, vh, ph, f1h, f2h, gh, uTop, uBot, vLef, vRig, au, bu, av, bv, infFlow, nu)

n = size(uh,1); h = 1 / n; 

nu = (infFlow==0)*nu + (infFlow>0)*h*infFlow/2;
rh1 = zeros(n,n+1); rh2 = zeros(n+1,n); rh3 = zeros(n,n);

%%
i = 2:n-1; j = 2:n;
rh1(i,j) = f1h(i,j) - (ph(i,j) - ph(i,j-1))/h + (nu*(uh(i-1,j) + uh(i+1,j) + uh(i,j-1) + uh(i,j+1) - 4*uh(i,j)) - ...
           0.5*h*au(i,j).*(uh(i,j+1) - uh(i,j-1)) - 0.5*h*bu(i,j).*(uh(i-1,j) - uh(i+1,j)))/h^2;
        
rh1(1,j) = f1h(1,j) - (ph(1,j) - ph(1,j-1))/h + (nu*(2*uTop(1,j) + uh(2,j) + uh(1,j-1) + uh(1,j+1) - 5*uh(1,j)) - ...
           0.5*h*au(1,j).*(uh(1,j+1) - uh(1,j-1)) - 0.5*h*bu(1,j).*(2*uTop(1,j) - uh(1,j) - uh(2,j)))/h^2;
        
rh1(n,j) = f1h(n,j) - (ph(n,j) - ph(n,j-1))/h + (nu*(uh(n-1,j) + 2*uBot(1,j) + uh(n,j-1) + uh(n,j+1) - 5*uh(n,j)) - ...
           0.5*h*au(n,j).*(uh(n,j+1) - uh(n,j-1)) - 0.5*h*bu(n,j).*(uh(n-1,j) + uh(n,j) - 2*uBot(1,j)))/h^2;
     
%%       
i = 2:n; j = 2:n-1;
rh2(i,j) = f2h(i,j) - (ph(i-1,j) - ph(i,j))/h + (nu*(vh(i-1,j) + vh(i+1,j) + vh(i,j-1) + vh(i,j+1) - 4*vh(i,j)) - ...
           0.5*h*av(i,j).*(vh(i,j+1) - vh(i,j-1)) - 0.5*h*bv(i,j).*(vh(i-1,j) - vh(i+1,j)))/h^2;
        
rh2(i,1) = f2h(i,1) - (ph(i-1,1) - ph(i,1))/h + (nu*(vh(i-1,1) + vh(i+1,1) + 2*vLef(i,1) + vh(i,2) - 5*vh(i,1)) - ...
           0.5*h*av(i,1).*(vh(i,1) + vh(i,2) - 2*vLef(i,1)) - 0.5*h*bv(i,1).*(vh(i-1,1) - vh(i+1,1)))/h^2;
        
rh2(i,n) = f2h(i,n) - (ph(i-1,n) - ph(i,n))/h + (nu*(vh(i-1,n) + vh(i+1,n) + vh(i,n-1) + 2*vRig(i,1) - 5*vh(i,n)) - ...
           0.5*h*av(i,n).*(2*vRig(i,1) - vh(i,n-1) - vh(i,n)) - 0.5*h*bv(i,n).*(vh(i-1,n) - vh(i+1,n)))/h^2;
    
%%
i = 1:n; j = 1:n;
rh3(i,j) = gh(i,j) + (uh(i,j+1) - uh(i,j))/h + (vh(i,j) - vh(i+1,j))/h;


end