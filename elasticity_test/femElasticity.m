function [err,soln] = femElasticity(mesh,pde,option,varargin)

%% Check input arguments
if isfield(mesh,'node') && isfield(mesh,'elem')
    node = mesh.node;
    elem = mesh.elem;
    h = mesh.h;
else
    warning(' mesh must contain node, elem and h')
end

if isfield(mesh,'bdFlag')
    bdFlag = mesh.bdFlag;
else
    warning(' mesh must need bdFlag'); 
end

if ~exist('option','var'), option = []; end

if ~exist('pde','var')
    warning(' solver need pde'); 
end

%% Parameters
nv = size(elem,2);          
dim = size(node,2);
option = femoption(option);
elemType = option.elemType;     

%% Initialize err
if option.errorflag  
    errL2 = 0; errH1 = 0; 
end

%% Finite Element Method        
switch elemType
    case 'P1'     % piecewise linear function P1 element
%         soln = ElasticityP1(node,elem,bdFlag,pde,option);
        soln = elasticity(node,elem,bdFlag,pde,option);
    case 'Q1'     % piecewise linear function P1 element
        soln = PoissonQ1(node,elem,bdFlag,pde,option);
    case 'CR'     % piecewise linear function CR element
        soln = PoissonCR(node,elem,bdFlag,pde,option);
    case 'P2'     % piecewise quadratic function
        [soln, eqn] = PoissonP2(node,elem,bdFlag,pde,option);
    case 'P3'     % piecewise cubic function
        soln = PoissonP3(node,elem,bdFlag,pde,option);
    case 'WG'     % weak Galerkin element
        soln = PoissonWG(node,elem,bdFlag,pde,option);
end

if option.plotflag % show mesh and solution
    figure(1);
    showresult(node,elem,soln.u);
end

%% Compute error
if option.errorflag  
    % Du
    if isfield(pde,'Du1') && isfield(pde,'Du2')
        if isfield(soln,'Du1') && ~isempty(soln.Du1) ...
           && isfield(soln,'Du2') && ~isempty(soln.Du2)
                % Du is in the output
            errH1u1 = getH1error(node,elem,pde.Du1,soln.Du1);
            errH1u2 = getH1error(node,elem,pde.Du2,soln.Du2);
        else
            errH1u1 = getH1error(node,elem,pde.Du1,soln.u1);          
            errH1u2 = getH1error(node,elem,pde.Du2,soln.u2);          
        end
        errH1   = sqrt(errH1u1^2 + errH1u2^2);
    end
    % L2
    if isfield(pde,'exactu1') && isfield(pde,'exactu2') 
        uh1 = soln.u1;
        uh2 = soln.u2;
        errL2u1 = getL2error(node,elem,pde.exactu1,uh1);
        errL2u2 = getL2error(node,elem,pde.exactu2,uh2);
        errL2 = sqrt(errL2u1^2+errL2u2^2);
%         % interpolation
%         switch elemType
%             case 'P1'
%                 uI = Lagrangeinterpolate(pde.exactu,node,elem);
%             case 'Q1'
%                 uI = Lagrangeinterpolate(pde.exactu,node,elem);
%             case 'CR'
%                 uI = Lagrangeinterpolate(pde.exactu,node,elem,'CR',eqn.edge);
%             case 'P2'
%                 uI = Lagrangeinterpolate(pde.exactu,node,elem,'P2',eqn.edge);
%             case 'P3'
%                 uI = Lagrangeinterpolate(pde.exactu,node,elem,'P3',eqn.edge);
%             case 'WG'
%                 uI = Lagrangeinterpolate(pde.exactu,node,elem,'WG',eqn.edge);
%         end
%         erruIuh = sqrt((uh-uI)'*eqn.A*(uh-uI));
%         errMax = max(abs(uh-uI));
    end

    N = length(soln.u);
    err = struct('h',h,'N',N,'H1',errH1,'L2',errL2 ); %,...
        % 'uIuhH1',erruIuh,'uIuhMax',errMax);
else
    err =[];
end         
            
