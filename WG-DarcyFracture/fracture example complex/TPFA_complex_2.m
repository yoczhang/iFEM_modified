function r = TPFA_complex_2(NGamma, lGamma, K, fracLength, ite)
% for the case ite is great than 1, ite = 2, 3, ...

co = lGamma.*K; alpha = 2*co./fracLength;
r = sparse(NGamma, NGamma); 

frac = cell(10,1);
frac{1} = (1:12*ite)'; % 1st is local num, 2nd is global num
frac{2} = 12*ite+(1:8*ite)';
frac{3} = 20*ite+(1:15*ite)';
frac{4} = 35*ite+(1:12*ite)';
frac{5} = 47*ite+(1:17*ite)';
frac{6} = 64*ite+(1:4*ite)';
frac{7} = 68*ite+(1:7*ite)';
frac{8} = 75*ite+(1:14*ite)';
frac{9} = 89*ite+(1:6*ite)';
frac{10} = 95*ite+(1:8*ite)';

%% Left
left = [frac{1}(1) frac{2}(1) frac{3}(1) frac{4}(1) frac{5}(1) ...
        frac{6}(1) frac{7}(1) frac{8}(1) frac{9}(1) frac{10}(1)];
for i = 1:length(left)
    j = left(i);
    r(j,j)   = r(j,j)   + alpha(j)*alpha(j+1)/(alpha(j)+alpha(j+1));
    r(j,j+1) = r(j,j+1) - alpha(j)*alpha(j+1)/(alpha(j)+alpha(j+1));
end
clear left

%% Right
right = [frac{1}(12*ite) frac{2}(8*ite) frac{3}(15*ite) frac{4}(12*ite) ...
         frac{7}(7*ite)  frac{8}(14*ite) frac{9}(6*ite) frac{10}(8*ite)];
for i = 1:length(right)
    j = right(i);
    r(j,j)   = r(j,j)   + alpha(j)*alpha(j-1)/(alpha(j)+alpha(j-1));
    r(j,j-1) = r(j,j-1) - alpha(j)*alpha(j-1)/(alpha(j)+alpha(j-1));
end
clear right

%% Interior
flag1 = setdiff(1:12*ite, [1, 12*ite, 7*ite, 7*ite+1])';
flag2 = setdiff(1:8*ite, [1, 8*ite, 4*ite, 4*ite+1])';
flag3 = setdiff(1:15*ite, [1, 15*ite])';
flag4 = setdiff(1:12*ite, [1, 12*ite, 2*ite, 2*ite+1])';
flag5 = setdiff(1:17*ite, [1, ite, ite+1, 14*ite, 14*ite+1, 17*ite])';
flag6 = setdiff(1:4*ite, [1, 4*ite])';
flag7 = setdiff(1:7*ite, [1, 7*ite, 6*ite, 6*ite+1])';
flag8 = setdiff(1:14*ite, [1, 14*ite, ite, ite+1, 10*ite, 10*ite+1])';
flag9 = setdiff(1:6*ite, [1, 6*ite])';
flag10 = setdiff(1:8*ite, [1, 8*ite, ite, ite+1, 7*ite, 7*ite+1])'; 
interior = [frac{1}(flag1); frac{2}(flag2); frac{3}(flag3); frac{4}(flag4);
            frac{5}(flag5); frac{6}(flag6); frac{7}(flag7); frac{8}(flag8);
            frac{9}(flag9); frac{10}(flag10)];
for i = 1:length(interior)
    j =  interior(i);
    r(j,j)   = r(j,j)   + alpha(j)*alpha(j-1)/(alpha(j)+alpha(j-1)) + alpha(j)*alpha(j+1)/(alpha(j)+alpha(j+1));
    r(j,j-1) = r(j,j-1) - alpha(j)*alpha(j-1)/(alpha(j)+alpha(j-1));
    r(j,j+1) = r(j,j+1) - alpha(j)*alpha(j+1)/(alpha(j)+alpha(j+1));
end
clear interior flag1 flag2 flag3 flag4 flag5 flag6 flag7 flag8 flag9 flag10

%% Intesection
angle = [5.119253873298953e-01     1.610655443491531e+00     2.123410727191778e+00
         7.600410772429374e-01     8.564696224717313e-01     8.520317395310949e-01];

inte = [frac{2}(4*ite)  frac{1}(7*ite)  frac{2}(4*ite+1) frac{1}(7*ite+1);
        frac{4}(2*ite)  frac{10}(ite+1) frac{4}(2*ite+1) frac{10}(ite);
        frac{8}(ite)  frac{10}(7*ite+1) frac{8}(ite+1) frac{10}(7*ite);
        frac{8}(10*ite) frac{5}(ite) frac{8}(10*ite+1) frac{5}(ite+1);
        frac{7}(6*ite)  frac{5}(14*ite) frac{7}(6*ite+1) frac{5}(14*ite+1)];
for i = 1:size(inte,1)
    a = angle(i);
%     angle2 = [0 a 0 pi-a;
%               a 0 pi-a 0;
%               0 pi-a 0 a;
%               pi-a 0 a 0];
    angle2 = [0 a pi pi-a;
              a 0 pi-a pi;
              pi pi-a 0 a;
              pi-a pi a 0];
    allTran = alpha(inte(i,1)) + alpha(inte(i,2)) + alpha(inte(i,3)) + alpha(inte(i,4));
    for k = 1:4
        temp = setdiff(1:4, k);
        r(inte(i,k),inte(i,k)) = r(inte(i,k),inte(i,k)) + alpha(inte(i,k))*alpha(inte(i,temp(1)))/allTran*cos(angle2(k,temp(1))/2) + ...
                                                          alpha(inte(i,k))*alpha(inte(i,temp(2)))/allTran*cos(angle2(k,temp(2))/2) + ...
                                                          alpha(inte(i,k))*alpha(inte(i,temp(3)))/allTran*cos(angle2(k,temp(3))/2);
        r(inte(i,k),inte(i,temp(1))) = r(inte(i,k),inte(i,temp(1))) - ...
                                       alpha(inte(i,k))*alpha(inte(i,temp(1)))/allTran*cos(angle2(k,temp(1))/2);
        r(inte(i,k),inte(i,temp(2))) = r(inte(i,k),inte(i,temp(2))) - ...
                                       alpha(inte(i,k))*alpha(inte(i,temp(2)))/allTran*cos(angle2(k,temp(2))/2);
        r(inte(i,k),inte(i,temp(3))) = r(inte(i,k),inte(i,temp(3))) - ...
                                       alpha(inte(i,k))*alpha(inte(i,temp(3)))/allTran*cos(angle2(k,temp(3))/2);
    end
end  
clear inte angle2 allTran temp

Tran = alpha(frac{6}(4*ite))*alpha(frac{5}(17*ite))/(alpha(frac{6}(4*ite))+alpha(frac{5}(17*ite))) * cos(angle(6)/2);
r(frac{6}(4*ite),frac{6}(4*ite)) = r(frac{6}(4*ite),frac{6}(4*ite)) + Tran; 
r(frac{6}(4*ite),frac{5}(17*ite)) = r(frac{6}(4*ite),frac{5}(17*ite)) - Tran;
r(frac{5}(17*ite),frac{6}(4*ite)) = r(frac{5}(17*ite),frac{6}(4*ite)) - Tran; 
r(frac{5}(17*ite),frac{5}(17*ite)) = r(frac{5}(17*ite),frac{5}(17*ite)) + Tran;
clear Tran

left = [frac{2}(4*ite) frac{1}(7*ite) frac{4}(2*ite) frac{10}(7*ite) ...
        frac{8}(10*ite) frac{7}(6*ite) frac{5}(14*ite) frac{5}(17*ite) frac{6}(4*ite)];
right = [frac{2}(4*ite+1) frac{1}(7*ite+1) frac{4}(2*ite+1) frac{10}(ite+1) ...
         frac{8}(ite+1) frac{8}(10*ite+1) frac{5}(ite+1) frac{5}(14*ite+1)];
for i = 1:length(left)
    j = left(i);
    r(j,j)   = r(j,j)   + alpha(j)*alpha(j-1)/(alpha(j)+alpha(j-1));
    r(j,j-1) = r(j,j-1) - alpha(j)*alpha(j-1)/(alpha(j)+alpha(j-1));
end
for i = 1:length(right)
    j = right(i);
    r(j,j)   = r(j,j)   + alpha(j)*alpha(j+1)/(alpha(j)+alpha(j+1));
    r(j,j+1) = r(j,j+1) - alpha(j)*alpha(j+1)/(alpha(j)+alpha(j+1));
end
clear left right
end