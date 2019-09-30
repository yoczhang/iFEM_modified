function FractureDarcyP1P1RT1_P1DG_2(eipsilon, sigma, ite)

%% Mesh
basis_type = 1;
[node1,elem1,elem2edge1,edge1,fixedEdge1,freeEdge1,node2,elem2,elem2edge2, ...
       edge2,fixedEdge2,freeEdge2,brother,Pb,Tb] = formMeshInformationFracturedDarcy_2(ite, basis_type);
NE1 = size(edge1,1); NT1 = size(elem1,1); 
NE2 = size(edge2,1); NT2 = size(elem2,1); 
n = 2^(ite+1); dof1D = 2*n; % n is the number of element in the horizontal line

[gcrt, gprt] = generate_Gauss_reference_triangle(12);
[gcrl, gprl] = generate_Gauss_reference_1D(8);

TE1 = [2*elem2edge1(:,1)-1 2*elem2edge1(:,1)   2*elem2edge1(:,2)-1 ...
       2*elem2edge1(:,2)   2*elem2edge1(:,3)-1 2*elem2edge1(:,3)];
TE2 = [2*elem2edge2(:,1)-1 2*elem2edge2(:,1)   2*elem2edge2(:,2)-1 ...
       2*elem2edge2(:,2)   2*elem2edge2(:,3)-1 2*elem2edge2(:,3)];

%% Parameters
psi = 0.75; lGamma = 0.01; eitaGamma = lGamma/eipsilon;

%% System
[A1, B1, C1] = assembleStiffnessMatrixForPrimalDarcyP1P1RT1(node1, elem1, elem2edge1, edge1, gcrt, gprt, TE1);
[A2, B2, C2] = assembleStiffnessMatrixForPrimalDarcyP1P1RT1(node2, elem2, elem2edge2, edge2, gcrt, gprt, TE2);

b1 = assembleLoadVectorForPrimalDarcyP1P1RT1(node1, elem1, @f1, gcrt, gprt);
b2 = assembleLoadVectorForPrimalDarcyP1P1RT1(node2, elem2, @f2, gcrt, gprt);

[D, b3] = P1DG_2(1, sigma, -1, 1, n, @pGamma, @fGamma, psi, eitaGamma, lGamma, eipsilon, gcrl, gprl);
[I1, I2, I3, I4, I5, I6, I7, I8] = assembleStiffnessMatrixForPrimalDarcyInterfaceTermP1P1RT1P1DG_2(NE1, ...
                                                            NE2, dof1D, brother, n, psi, eitaGamma, ...
                                                            Tb, basis_type, node1, edge1, node2, edge2, gcrl, gprl);
                                                        
temp1 = sparse(3*NT1,3*NT2); temp2 = sparse(3*NT1,2*NE2); temp3 = sparse(3*NT1,dof1D);
temp4 = sparse(2*NE1,3*NT2); temp5 = sparse(3*NT2,dof1D);
bigA = [A1      B1      temp1   temp2   temp3;
        B1'     C1+I1   temp4   I2      I7;
        temp1'  temp4'  A2      B2      temp5;
        temp2'  I3      B2'     C2+I4   I8;
        temp3'  I5      temp5'  I6      D];
bigF = [b1; zeros(2*NE1,1); b2; zeros(2*NE2,1); b3];

%% Solve
fixedEdge1 = setdiff(fixedEdge1, brother(:,1));
Qp1b = Qb_uForPrimalDarcyP1P1RT1(node1, edge1, fixedEdge1, @p1);
fixedEdge1 = [2*fixedEdge1-1 2*fixedEdge1];

fixedEdge2 = setdiff(fixedEdge2, brother(:,2));
Qp2b = Qb_uForPrimalDarcyP1P1RT1(node2, edge2, fixedEdge2, @p2);
fixedEdge2 = [2*fixedEdge2-1 2*fixedEdge2];

bigu = zeros(3*NT1+2*NE1+3*NT2+2*NE2+dof1D,1);
for i = 1:size(fixedEdge1,1)
    j1 = fixedEdge1(i,1); j2 = fixedEdge1(i,2);
    bigu(3*NT1+j1,1) = Qp1b(i,1); bigu(3*NT1+j2,1) = Qp1b(i,2);
end
for i = 1:size(fixedEdge2,1)
    j1 = fixedEdge2(i,1);j2 = fixedEdge2(i,2);
    bigu(3*NT1+2*NE1+3*NT2+j1,1) = Qp2b(i,1);
    bigu(3*NT1+2*NE1+3*NT2+j2,1) = Qp2b(i,2);
end

bigF = bigF - bigA*bigu;
fixed = [3*NT1+fixedEdge1(:); 3*NT1+2*NE1+3*NT2+fixedEdge2(:)];
free = setdiff(1:3*NT1+2*NE1+3*NT2+2*NE2+dof1D, fixed);
bigu(free) = bigA(free,free) \ bigF(free);

p10 = bigu(1:3*NT1); 
p20 = bigu(3*NT1+2*NE1+1:3*NT1+2*NE1+3*NT2); 

err1 = getErrorPrimalDarcyP1P1RT1(node1, elem1, gcrt, gprt, p10, @p1);
err2 = getErrorPrimalDarcyP1P1RT1(node2, elem2, gcrt, gprt, p20, @p2);
error = sqrt(err1^2 + err2^2);

pGammah = bigu(3*NT1+2*NE1+3*NT2+2*NE2+1:end);
[L2errPGamma, H1errPGamma] = getErrorFracturePartP1DG_2(pGammah, n, gcrl, gprl, @pGamma, @pGammax);
[error, L2errPGamma, H1errPGamma] 


%% Test for fracture
% Qp1b = Qb_uForPrimalDarcyP1P1RT1_allEdge(node1, edge1, (1:NE1)', @p1, gcrl, gprl);
% Qp2b = Qb_uForPrimalDarcyP1P1RT1_allEdge(node2, edge2, (1:NE2)', @p2, gcrl, gprl);
% b3 = b3 - I5*Qp1b - I6*Qp2b;
% pGammah = D \ b3;
% [L2errPGamma, H1errPGamma] = getErrorFracturePartP1DG_2(pGammah, n, gcrl, gprl, @pGamma, @pGammax);
% [L2errPGamma, H1errPGamma]
%% Data
    function s = p1(p)
        x = p(:,1); y = p(:,2);
        s = eipsilon*cos(x).*cosh(y) + (1-eipsilon)*cosh(lGamma/2)*cos(x);
    end
    function s = p2(p)
        x = p(:,1); y = p(:,2);
        s = eipsilon*cos(x).*cosh(y) + (1-eipsilon)*cosh(lGamma/2)*cos(x);
    end
    function s = pGamma(x)
        s = eipsilon*cos(x) + (1-eipsilon)*cosh(lGamma/2)*cos(x);
    end
    function s = pGammax(x)
        s = cosh(lGamma/2)*sin(x)*(eipsilon - 1) - eipsilon*sin(x);
    end
    function s = f1(p)
        x = p(:,1); y = p(:,2);
        s = (1-eipsilon)*cosh(lGamma/2)*cos(x);
    end
    function s = f2(p)
        x = p(:,1); y = p(:,2);
        s = (1-eipsilon)*cosh(lGamma/2)*cos(x);
    end
    function s = fGamma(x)
        s = eipsilon^2*cos(x) + eipsilon*(1-eipsilon)*cosh(lGamma/2)*cos(x);
    end
end