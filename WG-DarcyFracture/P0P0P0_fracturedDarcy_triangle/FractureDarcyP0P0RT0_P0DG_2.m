function [error0, errorb, H1p, H1pGamma, maxpGamma, dof]  = FractureDarcyP0P0RT0_P0DG_2(eipsilon, ite)

%% Mesh
basis_type = 1;
[node1,elem1,elem2edge1,edge1,fixedEdge1,freeEdge1,node2,elem2,elem2edge2, ...
       edge2,fixedEdge2,freeEdge2,brother,Pb,Tb] = formMeshInformationFracturedDarcy_2(ite, basis_type);
NE1 = size(edge1,1); NT1 = size(elem1,1); 
NE2 = size(edge2,1); NT2 = size(elem2,1); 
n = 2^ite; % n is the number of element in the horizontal line

[gcrt, gprt] = generate_Gauss_reference_triangle(9);
[gcrl, gprl] = generate_Gauss_reference_1D(8);

%% Parameters
psi = 0.75; lGamma = 0.01; eitaGamma = lGamma/eipsilon;

%% System
[A1, B1, C1, area1] = assembleStiffnessMatrixForPrimalDarcyP0P0RT0(node1, elem1, elem2edge1, NT1, NE1);
[A2, B2, C2, area2] = assembleStiffnessMatrixForPrimalDarcyP0P0RT0(node2, elem2, elem2edge2, NT2, NE2);

b1 = assembleLoadVectorForPrimalDarcyP0P0RT0(node1, elem1, NE1, @f1, gcrt, gprt);
b2 = assembleLoadVectorForPrimalDarcyP0P0RT0(node2, elem2, NE2, @f2, gcrt, gprt);

[D1, D2, b3] = P0DG_2(n, 1, lGamma, eipsilon, psi, eitaGamma, @fGamma, gcrl, gprl, pGamma(0), pGamma(1));
[I1, I2, I3, I4, I5, I6, I7, I8] = assembleStiffnessMatrixForPrimalDarcyInterfaceTermP0P0RT0P0DG_2(NE1, ...
                                                            NE2, brother, n, psi, eitaGamma);
                                                        
temp1 = sparse(NT1,NT2); temp2 = sparse(NT1,NE2); temp3 = sparse(NT1,n);
temp4 = sparse(NE1,NT2); temp5 = sparse(NT2,n);
bigA = [A1      B1      temp1   temp2   temp3;
        B1'     C1+I1   temp4   I2      I7;
        temp1'  temp4'  A2      B2      temp5;
        temp2'  I3      B2'     C2+I4   I8;
        temp3'  I5      temp5'  I6      D1+D2];
bigF = [b1; zeros(NE1,1); b2; zeros(NE2,1); b3];

%% Solve
[Qbp1, edgeLength1] = Qb_uForPrimalDarcyP0P0RT0(node1, edge1, (1:NE1)', @p1, gcrl, gprl);
[Qbp2, edgeLength2] = Qb_uForPrimalDarcyP0P0RT0(node2, edge2, (1:NE2)', @p2, gcrl, gprl);

fixedEdge1 = setdiff(fixedEdge1, brother(:,1));
fixedEdge2 = setdiff(fixedEdge2, brother(:,2));

dof = NT1+NE1+NT2+NE2+n;
bigu = zeros(NT1+NE1+NT2+NE2+n,1);
for i = 1:length(fixedEdge1)
    j = fixedEdge1(i); bigu(NT1+j,1) = Qbp1(j);
end
for i = 1:length(fixedEdge2)
    j = fixedEdge2(i); bigu(NT1+NE1+NT2+j,1) = Qbp2(j);
end

bigF = bigF - bigA*bigu;
fixed = [NT1+fixedEdge1; NT1+NE1+NT2+fixedEdge2];
free = setdiff(1:NT1+NE1+NT2+NE2+n, fixed);
% bigu(free) = bigA(free,free) \ bigF(free);
bigu(free) = amg(bigA(free,free), bigF(free));
% condest(bigA(free,free))

p10 = bigu(1:NT1); p20 = bigu(NT1+NE1+1:NT1+NE1+NT2); 
Q0p1 = Q0_uForPrimalDarcyP0P0RT0(node1, elem1, @p1, area1);
Q0p2 = Q0_uForPrimalDarcyP0P0RT0(node2, elem2, @p2, area2);
ph0 = [p10; p20]; Q0p = [Q0p1; Q0p2]; area = [area1; area2];
error0 = sqrt(sum(area.*((ph0 - Q0p).^2),1));
% figure; pdesurf(node1', elem1', p10'); shading flat; colorbar; hold on;
%         pdesurf(node2', elem2', p20'); shading flat; colorbar; hold off;

phb = [bigu(NT1+1:NT1+NE1); bigu(NT1+NE1+NT2+1:NT1+NE1+NT2+NE2)];
Qbp = [Qbp1; Qbp2]; edgeLength = [edgeLength1; edgeLength2];
hK = 2*sqrt(2)/n; errorb = sqrt(hK*sum(edgeLength.*((phb - Qbp).^2),1));

value1 = computeGradientError(node1, elem1, elem2edge1, area1, p10-Q0p1, phb(1:NE1)-Qbp1);
value2 = computeGradientError(node2, elem2, elem2edge2, area2, p20-Q0p2, phb(NE1+1:end)-Qbp2);
H1p = sqrt(value1^2 + value2^2);

pGammah = bigu(NT1+NE1+NT2+NE2+1:end);
h = 1/n; middle = (h/2: h: 1-h/2)';
pGammaI = pGamma(middle);

pGammae = pGammaI - pGammah;
H1pGamma = sqrt(pGammae'*D1*pGammae);
maxpGamma = max(abs(pGammae));

%% Data
    function s = p1(p)
        x = p(:,1); y = p(:,2);
%         s = eipsilon*cos(x).*cosh(y) + (1-eipsilon)*cosh(lGamma/2)*cos(x);
        s = (1-eipsilon)*exp(x+y) + lGamma*(1-eipsilon)/eipsilon*exp(x) + exp(y) + eitaGamma;
    end
    function s = p2(p)
        x = p(:,1); y = p(:,2);
%         s = eipsilon*cos(x).*cosh(y) + (1-eipsilon)*cosh(lGamma/2)*cos(x);
        s = (1-eipsilon)*exp(x+y) + exp(y);
    end
    function s = pGamma(x)
%         s = eipsilon*cos(x) + (1-eipsilon)*cosh(lGamma/2)*cos(x);
        s = (1-eipsilon)*exp(x) + lGamma*(1-eipsilon)/(2*eipsilon)*exp(x) + 1 + eitaGamma/2;
    end
    function s = pGammax(x)
%         s = cosh(lGamma/2)*sin(x)*(eipsilon - 1) - eipsilon*sin(x);
        s = -(exp(x)*(2*eipsilon + lGamma)*(eipsilon - 1))/(2*eipsilon);
    end
    function s = f1(p)
        x = p(:,1); y = p(:,2);
%         s = (1-eipsilon)*cosh(lGamma/2)*cos(x);
        s = 2*exp(x + y)*(eipsilon - 1) - exp(y) + (lGamma*exp(x)*(eipsilon - 1))/eipsilon;
    end
    function s = f2(p)
        x = p(:,1); y = p(:,2);
%         s = (1-eipsilon)*cosh(lGamma/2)*cos(x);
        s = 2*exp(x + y)*(eipsilon - 1) - exp(y);
    end
    function s = fGamma(x)
%         s = eipsilon^2*cos(x) + eipsilon*(1-eipsilon)*cosh(lGamma/2)*cos(x);
        s = (exp(x)*(2*eipsilon + lGamma)*(eipsilon - 1))/2;
    end
end