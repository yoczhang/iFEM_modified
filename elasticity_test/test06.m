%% Test elasticity equation
% Author:   Qi Li
% Modified: 2019/02/26

close all;
clear; clc;

addpath('./tool');

%% Parametere and data
para.mu = 1;
para.lambda = 1;
pde = elasticitydata06(para);

%% Space step
hh   =  (1/2).^(0:5)';  % test
maxIt = size(hh,1);

%% Type of basis function
option.elemType = 'P1';
% option.elemType = 'P2';
option.plotflag =0;
option.plotmesh = 0;

%% Initialization errors and orders
if isfield(pde,'exactu1') && isfield(pde,'exactu2') 
    option.errorflag = 1;    
    N = zeros(maxIt,1);
    H1error = zeros(maxIt,1);
    L2error = zeros(maxIt,1);
    H1orders = zeros(maxIt,1);
    L2orders = zeros(maxIt,1);
else
    option.errorflag = 0; 
end

for k=1:maxIt    
    h = hh(k);
   
%    domain = [0,1,0,1]; % domian 1
    domain = [-1,1,-1,1]; % domian 2
%    domain = [0,2,0,1]; % domian 3
%    domain = [-1,0,-1,0]; % domian 4      

    left  =  domain(1);
    right =  domain(2);
    bottom = domain(3);
    top =    domain(4);
    [node,elem] = squaremesh([left,right,bottom,top],h); % ok
    
    %% Dirichlet boundary. Result is good!
%     bdFlag = setboundary(node,elem,'Dirichlet');

    %% right Neumann boundary. Result is bad!
%     bdFlag = setboundary(node,elem,'Dirichlet','~(x==1)','Neumann','x==1');   
    bdFlag = setboundary(node,elem,'Dirichlet','all');   


    if 1 == option.plotmesh
        showmesh(node,elem);                                % plot mesh
        findelem(node,elem,'all','index','FaceColor','g');  % plot element indices
        findnode(node,'all','index','color','r');           % plot node indices
    end
    
    mesh=struct('node',node,'elem',elem,'bdFlag',bdFlag,'h',h);   
    
    %% Solve
    if 1 == option.errorflag 
        [err(k),solu] = femElasticity(mesh,pde,option);
        N(k) = err(k).N;
    else
        [~,solu] = femElasticity(mesh,pde,option);
    end    
    
    if k~=maxIt
%     [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
end

%% Convergence order
if option.errorflag     
    for i = 1:maxIt-1
        H1error(i) = err(i).H1;
        L2error(i) = err(i).L2;
        H1orders(i+1,1) = log(err(i).H1./err(i+1).H1)/log(hh(i)/hh(i+1));
        L2orders(i+1,1) = log(err(i).L2./err(i+1).L2)/log(hh(i)/hh(i+1));        
    end
    H1error(maxIt) = err(maxIt).H1;
    L2error(maxIt) = err(maxIt).L2;
    
    fprintf('\nRunning %s elemType=%s \n',pde.name,option.elemType);
    colname = {'#Dof','h','||u-u_h||_{H^1}','order','||u-u_h||_{L^2}','order'};
    disptable(colname,N,'%d',[ones(size(hh)),1./hh],'%d/%d', ...
              H1error, '%0.4e',H1orders, '%0.2f', ...
              L2error, '%0.4e',L2orders, '%0.2f');
end

%% Compute deformed mesh 
u = zeros(size(solu.u));
u(1:2:end) = solu.u1; u(2:2:end) = solu.u2;
[AvE,Eps3,Eps4,AvS,Sigma3,Sigma4] = ...
    avmatrix(node,elem,[],u,pde.lambda,pde.mu);

estimate = aposteriori(node,elem,[],AvE,Eps3,Eps4, ...
    AvS,Sigma3,Sigma4,u,pde.lambda,pde.mu)

%% Show results.
if size(elem,1)<=2500
    %% Show the mesh and boundary nodes.    
    subplot(2,1,1)
    showmesh(node,elem);
%     findnode(node,'all','index','color','r');
    
    %% Plot deformed mesh
    subplot(2,1,2)
    show(elem,[],node,AvS,u,pde.lambda,pde.mu,0.1);
end

% Dirichlet   good!
% Running elasticitydata01 elemType=P1 
%  #Dof      h    ||u-u_h||_{H^1} order ||u-u_h||_{L^2} order 
% 
%    162     1/4   2.4504e+00   0.00   3.4965e-01   0.00
%    578     1/8   1.2391e+00   0.98   1.0494e-01   1.74
%   2178    1/16   6.1795e-01   1.00   2.7711e-02   1.92
%   8450    1/32   3.0858e-01   1.00   7.0295e-03   1.98
%  33282    1/64   1.5423e-01   1.00   1.7639e-03   1.99
% 132098   1/128   7.7109e-02   1.00   4.4139e-04   2.00
% 526338   1/256   3.8553e-02   1.00   1.1037e-04   2.00


% Neumann      bad!
% Running elasticitydata06 elemType=P1 
%  #Dof     h   ||u-u_h||_{H^1} order ||u-u_h||_{L^2} order 
% 
%   18    1/1   6.7738e+00   0.00   8.8750e-01   0.00
%   50    1/2   4.3732e+00   0.63   7.9909e-01   0.15
%  162    1/4   2.7137e+00   0.69   3.6930e-01   1.11
%  578    1/8   1.8321e+00   0.57   2.4178e-01   0.61
% 2178   1/16   1.5342e+00   0.26   2.4092e-01   0.01

