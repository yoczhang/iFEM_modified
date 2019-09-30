function [L2err, H1err] = getErrorFracturePartP1DG_2(uh, n, gcrl,gprl, u, ux)

L2err = 0; H1err = 0; h = 2/n; left = -1; right = 1;
middle = left+h/2 : h : right-h/2; 
for i = 1:n
    uhLocal = uh(2*i-1:2*i);
    [gcll,gpll] = generate_Gauss_local_1D(gcrl,gprl,[left+(i-1)*h,0], [left+i*h,0], 0);
    
    for k = 1:length(gcll)
        L2err = L2err + gcll(k)*(u(gpll(k))  - basisDG1D(gpll(k),middle(i),0)*uhLocal)^2;
        H1err = H1err + gcll(k)*(ux(gpll(k)) - basisDG1D(gpll(k),middle(i),1)*uhLocal)^2;
    end
end
L2err = sqrt(L2err); H1err  = sqrt(H1err); 

    function s = basisDG1D(x, xc, order)
        if order==0
            s = [1+0*(x-xc), x-xc]; % xc is the barycenter
        elseif order==1
            s = [0+0*(x-xc), 1+0*(x-xc)];
        end
    end
end