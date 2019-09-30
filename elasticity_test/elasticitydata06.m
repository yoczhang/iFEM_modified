function pde = elasticitydata06(para)
%% ELASTICITYDATA data for elasticity problem
% 
% - mu \Delta u - (mu + lambda) grad(div u) = f    in \Omega
%                                         u = g_D  on \partial \Omega
%
% Created by Qi Li Wednesday, 27 Feb 2019.
%
%---------------------------------------------------------------------

%% Lame constants
if nargin == 0
    lambda = 1;
    mu = 1;
else
    if ~isstruct(para)
        exit('we need a struct data');
    end
    if ~isfield(para,'lambda') || isempty(para.lambda)
        lambda = 1;
    else
        lambda = para.lambda;
    end
    if ~isfield(para,'mu') || isempty(para.mu)
        mu = 1;
    else
        mu = para.mu;
    end
end

pde = struct('lambda',lambda, ...
             'mu',mu, ...
             'f', @f, ...
             'exactu1',@exactu1, ...
             'exactu2',@exactu2, ...
             'Du1',@Du1, ...
             'Du2',@Du2, ...
             'g_D',@g_D, ...
             'g_N',@g_N, ...
             'name','elasticitydata06');
%% subfunctions %%%%%%
    % load data (right hand side function)
    function z = f(p)
       z = mu*2*pi^2.*[exactu1(p),exactu2(p)];
    end
    % exact solution u1
    function z = exactu1(p)
       x = p(:,1); y = p(:,2);
       z = cos(pi*x).*cos(pi*y);    
    end
    % exact solution u2
    function z = exactu2(p)
       x = p(:,1); y = p(:,2);
       z = sin(pi*x).*sin(pi*y);    
    end
    % Derivative of the exact solution u1
    function uprime = Du1(p)
        x = p(:,1); y = p(:,2);
        uprime(:,1) = -pi*cos(pi*y).*sin(pi*x);
        uprime(:,2) = -pi*cos(pi*x).*sin(pi*y);
    end
    % Derivative of the exact solution u2
    function uprime = Du2(p)
        x = p(:,1); y = p(:,2);
        uprime(:,1) = pi*cos(pi*x).*sin(pi*y);
        uprime(:,2) = pi*cos(pi*y).*sin(pi*x);
    end
    % Dirichlet boundary condition
    function z = g_D(p)
        z = [exactu1(p),exactu2(p)]; 
    end

    % Neumann boundary condition
    function g = g_N(p)
        
        %    domain = [0,1,0,1]; % domian 1
        domain = [-1,1,-1,1]; % domian 2
        %    domain = [0,2,0,1]; % domian 3
        %    domain = [-1,0,-1,0]; % domian 4
        
        left  =  domain(1);
        right =  domain(2);
        bottom = domain(3);
        top =    domain(4);
        
        g = zeros(size(p,1),2);
        x = p(:,1); y = p(:,2);
        uprime1 = [ -2*mu*pi*cos(pi*y).*sin(pi*x),zeros(size(p,1),1)];
        uprime2 = [ zeros(size(p,1),1), 2*mu*pi*cos(pi*y).*sin(pi*x)];

        leftbd = (abs(x-left)<eps);     leftbd_n = [-1,0]';  % left 
        g(leftbd,:) = [uprime1(leftbd,:)*leftbd_n,uprime2(leftbd,:)*leftbd_n];
        rightbd = (abs(x-right)<eps);    rightbd_n = [1,0]'; % right
        g(rightbd,:) = [uprime1(rightbd,:)*rightbd_n,uprime2(rightbd,:)*rightbd_n];
        topbd = (abs(y-top)<eps);        topbd_n = [0,1]';   % top
        g(topbd,:) = [uprime1(topbd,:)*topbd_n,uprime2(topbd,:)*topbd_n];
        bottombd = (abs(y-bottom)<eps);  bottombd_n = [0,-1]'; % bottom
        g(bottombd,:) = [uprime1(bottombd,:)*bottombd_n,uprime2(bottombd,:)*bottombd_n];

    end


end