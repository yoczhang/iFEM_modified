function r = TPFA_complex(NGamma, lGamma, K, fracLength, ite)

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
left = [frac{1}(1) frac{2}(1) frac{3}(1) frac{4}(1) frac{6}(1) frac{7}(1) frac{9}(1)];
for i = 1:length(left)
    j = left(i);
    r(j,j)   = r(j,j)   + alpha(j)*alpha(j+1)/(alpha(j)+alpha(j+1));
    r(j,j+1) = r(j,j+1) - alpha(j)*alpha(j+1)/(alpha(j)+alpha(j+1));
end
clear left

%% Right
right = [frac{1}(12) frac{2}(8) frac{3}(15) frac{4}(12) frac{8}(14) frac{9}(6)];
for i = 1:length(right)
    j = right(i);
    r(j,j)   = r(j,j)   + alpha(j)*alpha(j-1)/(alpha(j)+alpha(j-1));
    r(j,j-1) = r(j,j-1) - alpha(j)*alpha(j-1)/(alpha(j)+alpha(j-1));
end
clear right

%% Interior
interior = [frac{1}([2:6,9:11]); frac{2}([2,3,6,7]); frac{3}(2:14); ...
            frac{4}(4:11); frac{5}([3:13,16]); frac{6}([2,3]); ...
            frac{7}(2:5); frac{8}([3:9,12,13]); frac{9}(2:5); frac{10}(3:6)];
for i = 1:length(interior)
    j =  interior(i);
    r(j,j)   = r(j,j)   + alpha(j)*alpha(j-1)/(alpha(j)+alpha(j-1)) + alpha(j)*alpha(j+1)/(alpha(j)+alpha(j+1));
    r(j,j-1) = r(j,j-1) - alpha(j)*alpha(j-1)/(alpha(j)+alpha(j-1));
    r(j,j+1) = r(j,j+1) - alpha(j)*alpha(j+1)/(alpha(j)+alpha(j+1));
end
clear interior

%% Intesection
angle = [5.119253873298953e-01     1.610655443491531e+00     2.123410727191778e+00
         7.600410772429374e-01     8.564696224717313e-01     8.520317395310949e-01];

inte = [frac{2}(4)  frac{1}(7)  frac{2}(5) frac{1}(8);
        frac{4}(2)  frac{10}(2) frac{4}(3) frac{10}(1);
        frac{8}(1)  frac{10}(8) frac{8}(2) frac{10}(7);
        frac{8}(10) frac{5}(1) frac{8}(11) frac{5}(2);
        frac{7}(6)  frac{5}(14) frac{7}(7) frac{5}(15)];
for i = 1:size(inte,1)
    a = angle(i);
    angle2 = [0 a 0 pi-a;
              a 0 pi-a 0;
              0 pi-a 0 a;
              pi-a 0 a 0];
%     angle2 = [0 a pi pi-a;
%               a 0 pi-a pi;
%               pi pi-a 0 a;
%               pi-a pi a 0];
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

Tran = alpha(frac{6}(4))*alpha(frac{5}(17))/(alpha(frac{6}(4))+alpha(frac{5}(17))) * cos(angle(6)/2);
r(frac{6}(4),frac{6}(4)) = r(frac{6}(4),frac{6}(4)) + Tran; 
r(frac{6}(4),frac{5}(17)) = r(frac{6}(4),frac{5}(17)) - Tran;
r(frac{5}(17),frac{6}(4)) = r(frac{5}(17),frac{6}(4)) - Tran; 
r(frac{5}(17),frac{5}(17)) = r(frac{5}(17),frac{5}(17)) + Tran;
clear Tran

left = [frac{2}(4) frac{1}(7) frac{4}(2) frac{10}(7) frac{8}(10) frac{7}(6) frac{5}(14) frac{5}(17) frac{6}(4)];
right = [frac{2}(5) frac{1}(8) frac{4}(3) frac{10}(2) frac{8}(2) frac{8}(11) frac{5}(2) frac{5}(15)];
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
%% compute intersection angle
% v1 = node(edge(brother(frac{2}(4),1),1),:) - node(edge(brother(frac{2}(4),1),2),:);
% v2 = node(edge(brother(frac{1}(7),1),1),:) - node(edge(brother(frac{1}(7),1),2),:);
% angle1 = acos(v1*v2'/(norm(v1)*norm(v2))); clear v1 v2
% 
% v1 = node(edge(brother(frac{10}(1),1),1),:) - node(edge(brother(frac{10}(1),1),2),:);
% v2 = node(edge(brother(frac{4}(2),1),1),:) - node(edge(brother(frac{4}(2),1),2),:);
% angle2 = acos(v1*v2'/(norm(v1)*norm(v2))); clear v1 v2
% 
% v1 = node(edge(brother(frac{10}(7),1),1),:) - node(edge(brother(frac{10}(7),1),2),:);
% v2 = node(edge(brother(frac{8}(1),1),1),:) - node(edge(brother(frac{8}(1),1),2),:);
% angle3 = acos(v1*v2'/(norm(v1)*norm(v2))); clear v1 v2
% 
% v1 = node(edge(brother(frac{8}(10),1),1),:) - node(edge(brother(frac{8}(10),1),2),:);
% v2 = node(edge(brother(frac{5}(1),1),1),:) - node(edge(brother(frac{5}(1),1),2),:);
% angle4 = acos(v1*v2'/(norm(v1)*norm(v2))); clear v1 v2
% 
% v1 = node(edge(brother(frac{5}(14),1),1),:) - node(edge(brother(frac{5}(14),1),2),:);
% v2 = node(edge(brother(frac{7}(6),1),1),:) - node(edge(brother(frac{7}(6),1),2),:);
% angle5 = acos(v1*v2'/(norm(v1)*norm(v2))); clear v1 v2
% 
% v1 = node(edge(brother(frac{6}(4),1),1),:) - node(edge(brother(frac{6}(4),1),2),:);
% v2 = node(edge(brother(frac{5}(17),1),1),:) - node(edge(brother(frac{5}(17),1),2),:);
% angle6 = acos(v1*v2'/(norm(v1)*norm(v2))); clear v1 v2
% angle = [angle1 angle2 angle3 angle4 angle5 angle6]; 
% angle
end