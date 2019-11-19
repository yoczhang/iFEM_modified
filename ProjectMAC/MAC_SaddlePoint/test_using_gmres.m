function test_using_gmres()

clear
clc
close all

n = 11;
b = afun(ones(n,1));
tol = 1e-12;  maxit = 15; 
x1 = gmres(@afun,b,10,tol,maxit, @mfun);


    function y = afun(x)
        nn = length(x);
        y = [0; x(1:nn-1)] + ...
              [((nn-1)/2:-1:0)'; (1:(nn-1)/2)'].*x + ...
              [x(2:nn); 0];
        %y = 0*y;
    end

    function y = mfun(r)
        nn = length(r);
        y = r ./ [((nn-1)/2:-1:1)'; 1; (1:(nn-1)/2)'];
    end

end