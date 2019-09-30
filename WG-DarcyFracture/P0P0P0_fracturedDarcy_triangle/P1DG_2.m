function [A, b] = P1DG_2(symmetric, sigma, left, right, n, u, f, psi, eitaGamma, lGamma, eipsilon, gcrl,gprl)
% symmetric = -1. SIPG; 1, NIPG; 0, IIPG;
% sigma is the penalty parameter;
% n is the number of interval, so there are n+1 nodes;
% solve the problem:  - u_xx = f, [left,right]; u = uD, left or right.

%% Mesh
dis = right - left; h = dis/n; dof = 2*n; 
middle = left+h/2 : h : right-h/2; % barycenters
K = lGamma*eipsilon;
% [gcrl,gprl] = generate_Gauss_reference_1D(8);

%% Matrix
% Stiffness matrix
e = zeros(dof,1); e(2:2:dof) = h; A = spdiags(e, 0, dof, dof);  

% 1st node
value = temp(left, middle(1), middle(1), -1, -1, 1);
A(1,1) = A(1,1) + value(1,1); A(2,1) = A(2,1) + value(1,2);
A(1,2) = A(1,2) + value(2,1); A(2,2) = A(2,2) + value(2,2);

% last node
value = temp(right, middle(n), middle(n), 1, 1, 1);
A(2*n-1,2*n-1) = A(2*n-1,2*n-1) + value(1,1);
A(2*n-1,2*n)   = A(2*n-1,2*n)   + value(2,1);
A(2*n,2*n-1)   = A(2*n,2*n-1)   + value(1,2);
A(2*n,2*n)     = A(2*n,2*n)     + value(2,2);

% interior nodes
for i = 2:n
    node = left+(i-1)*h; % coordinate of current node
    eleLef = i-1; eleRig = i;
    xcLef = middle(eleLef); xcRig = middle(eleRig);
    
    % if test left and trial left
    value = temp(node, xcLef, xcLef, 1, 1, 0.5);
    A(2*eleLef-1,2*eleLef-1) = A(2*eleLef-1,2*eleLef-1) + value(1,1);
    A(2*eleLef,  2*eleLef-1) = A(2*eleLef,  2*eleLef-1) + value(1,2);
    A(2*eleLef-1,2*eleLef)   = A(2*eleLef-1,2*eleLef)   + value(2,1);
    A(2*eleLef,  2*eleLef)   = A(2*eleLef,  2*eleLef)   + value(2,2);
    
    % if test left and trial right
    value = temp(node, xcRig, xcLef, -1, 1, 0.5);
    A(2*eleLef-1,2*eleRig-1) = A(2*eleLef-1,2*eleRig-1) + value(1,1);
    A(2*eleLef,  2*eleRig-1) = A(2*eleLef,  2*eleRig-1) + value(1,2);
    A(2*eleLef-1,2*eleRig)   = A(2*eleLef-1,2*eleRig)   + value(2,1);
    A(2*eleLef,  2*eleRig)   = A(2*eleLef,  2*eleRig)   + value(2,2);
    
    % if test right and trial left
    value = temp(node, xcLef, xcRig, 1, -1, 0.5);
    A(2*eleRig-1,2*eleLef-1) = A(2*eleRig-1,2*eleLef-1) + value(1,1);
    A(2*eleRig,  2*eleLef-1) = A(2*eleRig,  2*eleLef-1) + value(1,2);
    A(2*eleRig-1,2*eleLef)   = A(2*eleRig-1,2*eleLef)   + value(2,1);
    A(2*eleRig,  2*eleLef)   = A(2*eleRig,  2*eleLef)   + value(2,2);
    
    % if test right and trial right
    value = temp(node, xcRig, xcRig, -1, -1, 0.5);
    A(2*eleRig-1,2*eleRig-1) = A(2*eleRig-1,2*eleRig-1) + value(1,1);
    A(2*eleRig,  2*eleRig-1) = A(2*eleRig,  2*eleRig-1) + value(1,2);
    A(2*eleRig-1,2*eleRig)   = A(2*eleRig-1,2*eleRig)   + value(2,1);
    A(2*eleRig,  2*eleRig)   = A(2*eleRig,  2*eleRig)   + value(2,2);
end
    function s = temp(x, xcTrial, xcTest, flagTrial, flagTest, co)
        s = -K*co*basisDG1D(x, xcTrial, 1)'*basisDG1D(x, xcTest, 0)*flagTest + ...
            K*symmetric*co*flagTrial*basisDG1D(x, xcTrial, 0)'*basisDG1D(x, xcTest, 1) + ...
            sigma/h*flagTrial*basisDG1D(x, xcTrial, 0)'*flagTest*basisDG1D(x, xcTest, 0);
    end

% mass matrix
e = repmat([h, h^3/12], 1, n);
B = spdiags(e', 0, dof, dof);

A = A + 4/((2*psi-1)*eitaGamma)*B;

%% Load vector
b = zeros(dof,1);
for i = 1:n
    [gcll,gpll] = generate_Gauss_local_1D(gcrl,gprl,[left+(i-1)*h,0], [left+i*h,0], 0);
    value = zeros(1,2);
    for k = 1:length(gcll)
        value = value + lGamma*gcll(k)*f(gpll(k))*basisDG1D(gpll(k), middle(i), 0);
    end
    b(2*i-1) = value(1); b(2*i) = value(2);
end
b(1) = b(1) + sigma/h*u(left); 
b(2) = b(2) - K*symmetric*u(left) + sigma/h*u(left)*(-h/2);
b(2*n-1) = b(2*n-1) + sigma/h*u(right);
b(2*n)   = b(2*n)   + K*symmetric*u(right) + sigma/h*u(right)*h/2;


    function s = basisDG1D(x, xc, order)
        if order==0
            s = [1+0*(x-xc), x-xc]; % xc is the barycenter
        elseif order==1
            s = [0+0*(x-xc), 1+0*(x-xc)];
        end
    end
end