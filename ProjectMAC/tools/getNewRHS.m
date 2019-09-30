function [r4, r5, r6] = getNewRHS(r1, r2, r3, alpha_u, alpha_p, au, bu, av, bv, infFlow)
 
% alpha_u = (infFlow==0)*1 + (infFlow>0)*alpha_u;
% alpha_p = (infFlow==0)*1 + (infFlow>0)*alpha_p;
% nu = (infFlow==0)*0 + (infFlow>0)*0.5*infFlow/size(r1,1);

nu = 0.5*infFlow/size(r1,1);
[t4, t5, t6] = getRHS(nu/2, r1, r2, r3);
t4 = alpha_u*t4; t5 = alpha_u*t5; t6 = alpha_p*t6;

[t1, t2, t3] = getRHS(nu, r1, r2, r3);

r4 = t1 - t4; r5 = t2 - t5; r6 = t3 - t6;

%% Subfunction to implement A*x in matrix-free way
    function [b1, b2, b3] = getRHS(nu, e1, e2, e3)
        m = size(e1,1); H = 1 / m; 
        b1 = zeros(m,m+1); b2 = zeros(m+1,m); b3 = zeros(m,m);
        
        i = 2:m-1; j = 2:m;
        b1(i,j) = nu*(4*e1(i,j) - e1(i-1,j) - e1(i+1,j) - e1(i,j-1) - e1(i,j+1))/H^2 + (0.5*au(i,j).*(e1(i,j+1) - e1(i,j-1)) + ...
                   0.5*bu(i,j).*(e1(i-1,j) - e1(i+1,j)))/H + (e3(i,j) - e3(i,j-1))/H;
                
        b1(1,j) = nu*(5*e1(1,j) - e1(2,j) - e1(1,j-1) - e1(1,j+1))/H^2 + (0.5*au(1,j).*(e1(1,j+1) - e1(1,j-1)) + ...
                  0.5*bu(1,j).*(-e1(1,j) - e1(2,j)))/H + (e3(1,j) - e3(1,j-1))/H;
                 
        b1(m,j) = nu*(5*e1(m,j) - e1(m-1,j) - e1(m,j-1) - e1(m,j+1))/H^2 + (0.5*au(m,j).*(e1(m,j+1) - e1(m,j-1)) + ...
                   0.5*bu(m,j).*(e1(m-1,j) + e1(m,j)))/H + (e3(m,j) - e3(m,j-1))/H;
         

        i = 2:m; j = 2:m-1;
        b2(i,j) = nu*(4*e2(i,j) - e2(i-1,j) - e2(i+1,j) - e2(i,j-1) - e2(i,j+1))/H^2 + (0.5*av(i,j).*(e2(i,j+1) - e2(i,j-1)) + ...
                   0.5*bv(i,j).*(e2(i-1,j) - e2(i+1,j)))/H + (e3(i-1,j) - e3(i,j))/H;
                
        b2(i,1) = nu*(5*e2(i,1) - e2(i-1,1) - e2(i+1,1) - e2(i,2))/H^2 + (0.5*av(i,1).*(e2(i,1) + e2(i,2)) + ...
                   0.5*bv(i,1).*(e2(i-1,1) - e2(i+1,1)))/H + (e3(i-1,1) - e3(i,1))/H;
                
        b2(i,m) = nu*(5*e2(i,m) - e2(i-1,m) - e2(i+1,m) - e2(i,m-1))/H^2 + (0.5*av(i,m).*(-e2(i,m-1) - e2(i,m)) + ...
                   0.5*bv(i,m).*(e2(i-1,m) - e2(i+1,m)))/H + (e3(i-1,m) - e3(i,m))/H;
     
        i = 1:m; j = 1:m;
        b3(i,j) = - (e1(i,j+1) - e1(i,j))/H - (e2(i,j) - e2(i+1,j))/H;
    end

end