function [error, L2errPGamma, H1errPGamma] = FractureDarcyP0P0RT0_P1CG_2(eipsilon, ite)

%% Mesh
basis_type = 1;
[node1,elem1,elem2edge1,edge1,fixedEdge1,freeEdge1,node2,elem2,elem2edge2, ...
       edge2,fixedEdge2,freeEdge2,brother,Pb,Tb] = formMeshInformationFracturedDarcy_2(ite, basis_type);
NE1 = size(edge1,1); NT1 = size(elem1,1); 
NE2 = size(edge2,1); NT2 = size(elem2,1); 
n = 2^(ite+1); dof1D = length(Pb); % n is the number of element in the horizontal line

[gcrt, gprt] = generate_Gauss_reference_triangle(9);
[gcrl, gprl] = generate_Gauss_reference_1D(4);

%% Parameters
psi = 0.75; lGamma = 0.01; eitaGamma = lGamma/eipsilon;

%% System
[A1, B1, C1] = assembleStiffnessMatrixForPrimalDarcyP0P0RT0(node1, elem1, elem2edge1, NT1, NE1);
[A2, B2, C2] = assembleStiffnessMatrixForPrimalDarcyP0P0RT0(node2, elem2, elem2edge2, NT2, NE2);

b1 = assembleLoadVectorForPrimalDarcyP0P0RT0(node1, elem1, NE1, @f1, gcrt, gprt);
b2 = assembleLoadVectorForPrimalDarcyP0P0RT0(node2, elem2, NE2, @f2, gcrt, gprt);

[D, b3] = Poisson1D_2(Pb, Tb, basis_type, n, gcrl, gprl, lGamma, eipsilon, psi, eitaGamma, @fGamma);
[I1, I2, I3, I4, I5, I6, I7, I8] = assembleStiffnessMatrixForPrimalDarcyInterfaceTermP0P0RT0_2(NE1, ...
                                                            NE2, dof1D, brother, n, psi, eitaGamma, ...
                                                            Tb, basis_type);
                                                        
temp1 = sparse(NT1,NT2); temp2 = sparse(NT1,NE2); temp3 = sparse(NT1,dof1D);
temp4 = sparse(NE1,NT2); temp5 = sparse(NT2,dof1D);
bigA = [A1      B1      temp1   temp2   temp3;
        B1'     C1+I1   temp4   I2      I7;
        temp1'  temp4'  A2      B2      temp5;
        temp2'  I3      B2'     C2+I4   I8;
        temp3'  I5      temp5'  I6      D];
bigF = [b1; zeros(NE1,1); b2; zeros(NE2,1); b3];

%% Solve
Qbp1 = Qb_uForPrimalDarcyP0P0RT0(node1, edge1, (1:NE1)', @p1, gcrl, gprl);
Qbp2 = Qb_uForPrimalDarcyP0P0RT0(node2, edge2, (1:NE2)', @p2, gcrl, gprl);

fixedEdge1 = setdiff(fixedEdge1, brother(:,1));
fixedEdge2 = setdiff(fixedEdge2, brother(:,2));

bigu = zeros(NT1+NE1+NT2+NE2+dof1D,1);
for i = 1:length(fixedEdge1)
    j = fixedEdge1(i); bigu(NT1+j,1) = Qbp1(j);
end
for i = 1:length(fixedEdge2)
    j = fixedEdge2(i); bigu(NT1+NE1+NT2+j,1) = Qbp2(j);
end
bigu(NT1+NE1+NT2+NE2+1,1) = pGamma(Pb(1));
bigu(NT1+NE1+NT2+NE2+dof1D,1) = pGamma(Pb(dof1D));

bigF = bigF - bigA*bigu;
fixed = [NT1+fixedEdge1; NT1+NE1+NT2+fixedEdge2; NT1+NE1+NT2+NE2+1; NT1+NE1+NT2+NE2+dof1D];
free = setdiff(1:NT1+NE1+NT2+NE2+dof1D, fixed);
% bigu(free) = bigA(free,free) \ bigF(free); 
bigu(free) = amg(bigA(free,free), bigF(free));


% p1b = bigu(NT1+1:NT1+NE1);
% edgeDiff = p1b - Qbp1;
% abs(edgeDiff(brother(:,1)))
% abs(p1b - Qbp1)
% % pGammaI = pGamma(Pb');
% % bigu(NT1+NE1+NT2+NE2+1:NT1+NE1+NT2+NE2+dof1D) = pGammaI;
% % bigF = bigF - bigA*bigu;
% % fixed = [NT1+fixedEdge1; NT1+NE1+NT2+fixedEdge2; (NT1+NE1+NT2+NE2+1:NT1+NE1+NT2+NE2+dof1D)'];
% % free = setdiff(1:NT1+NE1+NT2+NE2+dof1D, fixed);
% % bigu(free) = bigA(free,free) \ bigF(free); 


p10 = bigu(1:NT1); 
% p1I = Q0_uForPrimalDarcyP0P0RT0(node1, elem1, gcrt, gprt, @p1);
p20 = bigu(NT1+NE1+1:NT1+NE1+NT2); 
% p2I = Q0_uForPrimalDarcyP0P0RT0(node2, elem2, gcrt, gprt, @p2);

% figure, pdesurf(node1', elem1', (p10-p1I)'); shading flat; colorbar; hold on;
% pdesurf(node2', elem2', (p20-p2I)'); shading flat; colorbar; hold off; view(2);
% figure, pdesurf(node1', elem1', p10'); shading flat; colorbar; hold on;
% pdesurf(node2', elem2', p20'); shading flat; colorbar; hold off;

err1 = getErrorPrimalDarcyP0P0RT0(node1, elem1, gcrt, gprt, p10, @p1);
err2 = getErrorPrimalDarcyP0P0RT0(node2, elem2, gcrt, gprt, p20, @p2);
error = sqrt(err1^2 + err2^2);

pGammah = bigu(NT1+NE1+NT2+NE2+1:end);
% pGammah = zeros(dof1D,1);
% rightTerm = b3 - I5*Qbp1 - I6*Qbp2;
% pGammah(1) = pGamma(Pb(1)); pGammah(dof1D) = pGamma(Pb(dof1D));
% rightTerm = rightTerm - D*pGammah;
% free = 2:dof1D-1;
% pGammah(free) = D(free,free) \ rightTerm(free);

[L2errPGamma, H1errPGamma] = getErrorFracturePart_2(pGammah, @pGamma, @pGammax, gcrl, gprl, Tb, basis_type);
% temp1 = (sin(2)*(125000*sinh(2) - (549212923164140125*sinh(1))/549755813888 + 28976035995029013217976622760641/38685626227668133590597632) + (atan(sin(1/2)/cos(1/2))*(4835703278458516698824704000000*sinh(2) - 38647423689204411217286791168000*sinh(1) + 28976035995029013217976622760641))/9671406556917033397649408)^(1/2);
% temp2 = (sin(2)*(125000*sinh(2) - (549212923164140125*sinh(1))/549755813888 + 28976035995029013217976622760641/38685626227668133590597632) + (atan(sin(1/2)/cos(1/2))*(4835703278458516698824704000000*sinh(2) - 38647423689204411217286791168000*sinh(1) + 28976035995029013217976622760641))/9671406556917033397649408)^(1/2);


% pGammaI = pGamma(Pb');
% figure, plot(Pb', pGammaI);
% figure, plot(Pb', pGammah);

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