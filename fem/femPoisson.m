function [err,time,solver,eqn] = femPoisson(mesh,pde,option,varargin)
%% FEMPOISSON solve Poisson equation by various finite element methods
%
%   FEMPOISSON computes approximations to the Poisson equation on a
%   sequence of meshes obtained by uniform refinement of a input mesh.
% 
% See also Poisson, crack, Lshape
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

%% Check input arguments
if isfield(mesh,'node') && isfield(mesh,'elem')
    node = mesh.node;
    elem = double(mesh.elem);
else
    [node,elem] = squaremesh([0,1,0,1],0.125);  % default mesh is a square
end
if isfield(mesh,'bdFlag')
    bdFlag = mesh.bdFlag;
else
    bdFlag = setboundary(node,elem,'Dirichlet'); 
end
if ~exist('option','var'), option = []; end
if ~exist('pde','var')
    pde = sincosdata;                          % default data
end

%% Parameters
nv = size(elem,2);          
dim = size(node,2);
option = femoption(option);
elemType = option.elemType;     
refType = option.refType;
maxIt = option.maxIt;       
maxN = option.maxN;         
L0 = option.L0;

%% Generate an initial mesh by refine L0 times
for k = 1:L0
    if strcmp(refType,'red')
        if nv == 4
            [node,elem,bdFlag] = uniformrefinequad(node,elem,bdFlag);
        else
            [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
        end        
    elseif strcmp(refType,'bisect')
        [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    end
end

%% Initialize err
errL2 = zeros(maxIt,1);   errH1 = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1); errMax = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1); h = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % refine mesh
    t = cputime;
%     if strcmp(refType,'red')
%         if nv == 4
%             [node,elem,bdFlag] = uniformrefinequad(node,elem,bdFlag);
%         else
%             [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%         end        
%     elseif strcmp(refType,'bisect')
%         [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
%     end
    meshTime(k) = cputime - t;    
    % solve the equation
    switch elemType
        case 'P1_periodic'     % piecewise linear function P1 element with periodic BC
            [soln,eqn,info] = PoissonPeriodic(node,elem,bdFlag,pde,option);
        case 'P1'     % piecewise linear function P1 element
            [soln,eqn,info] = Poisson(node,elem,bdFlag,pde,option);
        case 'Q1'     % piecewise linear function P1 element
            [soln,eqn,info] = PoissonQ1(node,elem,bdFlag,pde,option);
        case 'CR'     % piecewise linear function CR element
            [soln,eqn,info] = PoissonCR(node,elem,bdFlag,pde,option);
        case 'P2'     % piecewise quadratic function
            [soln,eqn,info] = PoissonP2(node,elem,bdFlag,pde,option);
        case 'P3'     % piecewise cubic function
            [soln,eqn,info] = PoissonP3(node,elem,bdFlag,pde,option);
        case 'WG'     % weak Galerkin element
            [soln,eqn,info] = PoissonWG(node,elem,bdFlag,pde,option);            
    end
    % compute error
    t = cputime;
    if isfield(pde,'Du') 
        if isfield(soln,'Du') && ~isempty(soln.Du) % Du is in the output
            errH1(k) = getH1error(node,elem,pde.Du,soln.Du);
        elseif (dim == 2) && (nv == 4) % bilinear element
            errH1(k) = getH1errorQ1(node,elem,pde.Du,soln.u);   
        else
            errH1(k) = getH1error(node,elem,pde.Du,soln.u);            
        end
    end
    if isfield(pde,'exactu')
        uh = soln.u;
        if (dim == 2) && (nv == 4) % bilinear element
            errL2(k) = getL2errorQ1(node,elem,pde.exactu,uh);
        else
            errL2(k) = getL2error(node,elem,pde.exactu,uh);
        end
        % interpolation
        switch elemType
            case 'P1'
                uI = Lagrangeinterpolate(pde.exactu,node,elem);
            case 'Q1'
                uI = Lagrangeinterpolate(pde.exactu,node,elem);
            case 'CR'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'CR',eqn.edge);
            case 'P2'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'P2',eqn.edge);
            case 'P3'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'P3',eqn.edge);
            case 'WG'
                uI = Lagrangeinterpolate(pde.exactu,node,elem,'WG',eqn.edge);
        end
        erruIuh(k) = sqrt((uh-uI)'*eqn.A*(uh-uI));
        errMax(k) = max(abs(uh-uI));
    end
    errTime(k) = cputime - t;
    % record time
    solverTime(k) = info.solverTime;
    assembleTime(k) = info.assembleTime;
    if option.printlevel>1
        fprintf('Time to compute the error %4.2g s \n H1 err %4.2g    L2err %4.2g \n',...
                 errTime(k), errH1(k), errL2(k));    
    end
    % record solver information
    itStep(k) = info.itStep;
    stopErr(k) = info.stopErr;
    flag(k) = info.flag;
    % plot 
    N(k) = length(soln.u);
    h(k) = 1./(sqrt(size(node,1))-1);
    if  strcmp(elemType,'WG') % modify size for WG
        if ~isfield(option,'reducesystem') || (option.reducesystem == 1)
            N(k) = N(k) - size(elem,1); % reduced system
        end    
    end            
    if option.plotflag % show mesh and solution
       figure(1);  
       showresult(node,elem,soln.u);  
    end
    if N(k) > maxN
        break;
    end
end

%% Plot convergence rates
if option.rateflag
    figure;
    set(gcf,'Units','normal'); 
    set(gcf,'Position',[0.25,0.25,0.55,0.4]);
    subplot(1,2,1)
    showrateh2(h(1:k),errH1(1:k),1,'-*','|| Du-Du_h||',...
               h(1:k),errL2(1:k),1,'k-+','|| u-u_h||');
    subplot(1,2,2)
    showrateh2(h(1:k),erruIuh(1:k),1,'m-+','|| Du_I-Du_h||',...
               h(1:k),errMax(1:k),1,'r-*','|| u_I-u_h||_{\infty}');
end

%% Output
err = struct('h',h(1:k),'N',N(1:k),'H1',errH1(1:k),'L2',errL2(1:k),...
             'uIuhH1',erruIuh(1:k),'uIuhMax',errMax(1:k));         
time = struct('N',N(1:k),'err',errTime(1:k),'solver',solverTime(1:k), ...
              'assemble',assembleTime(1:k),'mesh',meshTime(1:k));
solver = struct('N',N(1:k),'itStep',itStep(1:k),'time',solverTime(1:k),...
                'stopErr',stopErr(1:k),'flag',flag(1:k));
            
%% Display error and time
if option.dispflag
    fprintf('\n');
%     display('Table: Error')
    colname = {'#Dof','h','||u-u_h||','||Du-Du_h||','||DuI-Du_h||','||uI-u_h||_{max}'};
    disptable(colname,err.N,[],err.h,'%0.3e',err.L2,'%0.5e',err.H1,'%0.5e',...
                      err.uIuhH1,'%0.5e',err.uIuhMax,'%0.5e');

%     display('Table: CPU time')
    colname = {'#Dof','Assemble','Solve','Error','Mesh'};
    disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e',...
                      time.err,'%0.2e',time.mesh,'%0.2e');
    fprintf('\n');
end