function test_using_gmres()

n = 21;
b = afun(ones(n,1));
tol = 1e-12;  maxit = 15; 
x1 = gmres(@afun,b,10,tol,maxit);


    function y = afun(x)
        y = [0; x(1:n-1)] + ...
              [((n-1)/2:-1:0)'; (1:(n-1)/2)'].*x + ...
              [x(2:n); 0];
        y = 0*y;
    end

    function y = mfun(r)
        y = r ./ [((n-1)/2:-1:1)'; 1; (1:(n-1)/2)'];
    end

end