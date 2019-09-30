function r = TPFA_complex_refineMeshStep(NGamma, lGamma, K, fracLength, refineMeshStep)
% for any refineMeshStep, 0, 1, 2, 3, ...

co = lGamma.*K; a = 2*co./fracLength; % a is alpha, here for simplicity
r = sparse(NGamma, NGamma); 
     
frac = [12; 8; 15; 12; 17; 4; 7; 14; 6; 8]; frac = cumsum(frac) * 2^refineMeshStep;

%% Fracture 1
left = 1; i11 = 7*2^refineMeshStep; i12 = i11+1; righ = frac(1);
others = setdiff((1:frac(1))', [left; i11; i12; righ]);
temp1 = a(1)*a(2)/(a(1)+a(2)); r(1,1) = r(1,1) + temp1; r(1,2) = r(1,2) - temp1;
temp2 = a(i11)*a(i11-1)/(a(i11)+a(i11-1)); r(i11,i11) = r(i11,i11) + temp2; r(i11,i11-1) = r(i11,i11-1) - temp2;
temp3 = a(i12)*a(i12+1)/(a(i12)+a(i12+1)); r(i12,i12) = r(i12,i12) + temp3; r(i12,i12+1) = r(i12,i12+1) - temp3;
temp4 = a(righ)*a(righ-1)/(a(righ)+a(righ-1)); r(righ,righ) = r(righ,righ) + temp4; r(righ,righ-1) = r(righ,righ-1) - temp4;
for i = 1:length(others)
    j = others(i); temp5 = a(j)*a(j-1)/(a(j)+a(j-1)); temp6 = a(j)*a(j+1)/(a(j)+a(j+1));
    r(j,j) = r(j,j) + temp5 + temp6; r(j,j-1) = r(j,j-1) - temp5; r(j,j+1) = r(j,j+1) - temp6;
end
clear left righ others temp1 temp2 temp3 temp4 temp5 temp6

%% Fracture 2
left = frac(1)+1; i21 = frac(1)+4*2^refineMeshStep; i22 = i21 + 1; righ = frac(2);
others = setdiff((frac(1)+1:frac(2))', [left; i21; i22; righ]);
temp1 = a(left)*a(left+1)/(a(left)+a(left+1)); r(left,left) = r(left,left) + temp1; r(left,left+1) = r(left,left+1) - temp1;
temp2 = a(i21)*a(i21-1)/(a(i21)+a(i21-1)); r(i21,i21) = r(i21,i21) + temp2; r(i21,i21-1) = r(i21,i21-1) - temp2;
temp3 = a(i22)*a(i22+1)/(a(i22)+a(i22+1)); r(i22,i22) = r(i22,i22) + temp3; r(i22,i22+1) = r(i22,i22+1) - temp3;
temp4 = a(righ)*a(righ-1)/(a(righ)+a(righ-1)); r(righ,righ) = r(righ,righ) + temp4; r(righ,righ-1) = r(righ,righ-1) - temp4;
for i = 1:length(others)
    j = others(i); temp5 = a(j)*a(j-1)/(a(j)+a(j-1)); temp6 = a(j)*a(j+1)/(a(j)+a(j+1));
    r(j,j) = r(j,j) + temp5 + temp6; r(j,j-1) = r(j,j-1) - temp5; r(j,j+1) = r(j,j+1) - temp6;
end
clear left righ others temp1 temp2 temp3 temp4 temp5 temp6

%% Fracture 3
left = frac(2)+1; righ = frac(3); others = setdiff((frac(2)+1:frac(3))', [left; righ]);
temp1 = a(left)*a(left+1)/(a(left)+a(left+1)); r(left,left) = r(left,left) + temp1; r(left,left+1) = r(left,left+1) - temp1;
temp4 = a(righ)*a(righ-1)/(a(righ)+a(righ-1)); r(righ,righ) = r(righ,righ) + temp4; r(righ,righ-1) = r(righ,righ-1) - temp4;
for i = 1:length(others)
    j = others(i); temp5 = a(j)*a(j-1)/(a(j)+a(j-1)); temp6 = a(j)*a(j+1)/(a(j)+a(j+1));
    r(j,j) = r(j,j) + temp5 + temp6; r(j,j-1) = r(j,j-1) - temp5; r(j,j+1) = r(j,j+1) - temp6;
end
clear left righ others temp1 temp4 temp5 temp6

%% Fracture 9
left = frac(8)+1; righ = frac(9); others = setdiff((frac(8)+1:frac(9))', [left; righ]);
temp1 = a(left)*a(left+1)/(a(left)+a(left+1)); r(left,left) = r(left,left) + temp1; r(left,left+1) = r(left,left+1) - temp1;
temp4 = a(righ)*a(righ-1)/(a(righ)+a(righ-1)); r(righ,righ) = r(righ,righ) + temp4; r(righ,righ-1) = r(righ,righ-1) - temp4;
for i = 1:length(others)
    j = others(i); temp5 = a(j)*a(j-1)/(a(j)+a(j-1)); temp6 = a(j)*a(j+1)/(a(j)+a(j+1));
    r(j,j) = r(j,j) + temp5 + temp6; r(j,j-1) = r(j,j-1) - temp5; r(j,j+1) = r(j,j+1) - temp6;
end
clear left righ others temp1 temp4 temp5 temp6

%% Fracture 6
left = frac(5)+1; i62 = frac(6); others = setdiff((frac(5)+1:frac(6))', [left; i62]);
temp1 = a(left)*a(left+1)/(a(left)+a(left+1)); r(left,left) = r(left,left) + temp1; r(left,left+1) = r(left,left+1) - temp1;
temp4 = a(i62)*a(i62-1)/(a(i62)+a(i62-1)); r(i62,i62) = r(i62,i62) + temp4; r(i62,i62-1) = r(i62,i62-1) - temp4;
for i = 1:length(others)
    j = others(i); temp5 = a(j)*a(j-1)/(a(j)+a(j-1)); temp6 = a(j)*a(j+1)/(a(j)+a(j+1));
    r(j,j) = r(j,j) + temp5 + temp6; r(j,j-1) = r(j,j-1) - temp5; r(j,j+1) = r(j,j+1) - temp6;
end
clear left others temp1 temp4 temp5 temp6

%% Fracture 7
if refineMeshStep==0
    left = frac(6)+1; i72 = frac(7); i71 = i72-1; others = frac(6)+[2;3;4;5];
    temp1 = a(left)*a(left+1)/(a(left)+a(left+1)); r(left,left) = r(left,left) + temp1; r(left,left+1) = r(left,left+1) - temp1;
    temp2 = a(i71)*a(i71-1)/(a(i71)+a(i71-1)); r(i71,i71) = r(i71,i71) + temp2; r(i71,i71-1) = r(i71,i71-1) - temp2;
    for i = 1:length(others)
        j = others(i); temp5 = a(j)*a(j-1)/(a(j)+a(j-1)); temp6 = a(j)*a(j+1)/(a(j)+a(j+1));
        r(j,j) = r(j,j) + temp5 + temp6; r(j,j-1) = r(j,j-1) - temp5; r(j,j+1) = r(j,j+1) - temp6;
    end
    clear left others temp1 temp2 temp5 temp6
else
    left = frac(6)+1; i71 = frac(6)+6*2^refineMeshStep; i72 = i71+1; righ = frac(7); others = setdiff((frac(6)+1:frac(7))', [left; i71; i72; righ]);
    temp1 = a(left)*a(left+1)/(a(left)+a(left+1)); r(left,left) = r(left,left) + temp1; r(left,left+1) = r(left,left+1) - temp1;
    temp4 = a(righ)*a(righ-1)/(a(righ)+a(righ-1)); r(righ,righ) = r(righ,righ) + temp4; r(righ,righ-1) = r(righ,righ-1) - temp4;
    temp2 = a(i71)*a(i71-1)/(a(i71)+a(i71-1)); r(i71,i71) = r(i71,i71) + temp2; r(i71,i71-1) = r(i71,i71-1) - temp2;
    temp3 = a(i72)*a(i72+1)/(a(i72)+a(i72+1)); r(i72,i72) = r(i72,i72) + temp3; r(i72,i72+1) = r(i72,i72+1) - temp3;
    for i = 1:length(others)
        j = others(i); temp5 = a(j)*a(j-1)/(a(j)+a(j-1)); temp6 = a(j)*a(j+1)/(a(j)+a(j+1));
        r(j,j) = r(j,j) + temp5 + temp6; r(j,j-1) = r(j,j-1) - temp5; r(j,j+1) = r(j,j+1) - temp6;
    end
    clear left righ others temp1 temp2 temp3 temp4 temp5 temp6
end

%% Fracture 5
if refineMeshStep==0
    i51 = frac(4)+1; i52 = i51 + 1; left = []; 
else
    i51 = frac(4)+2^refineMeshStep; i52 = i51 + 1; left = frac(4) + 1;
    temp2 = a(i51)*a(i51-1)/(a(i51)+a(i51-1)); r(i51,i51) = r(i51,i51) + temp2; r(i51,i51-1) = r(i51,i51-1) - temp2; 
    temp1 = a(left)*a(left+1)/(a(left)+a(left+1)); r(left,left) = r(left,left) + temp1; r(left,left+1) = r(left,left+1) - temp1;
    clear temp1 temp2
end
i53 = frac(4) + 14*2^refineMeshStep; i54 = i53 + 1; i55 = frac(5);
others = setdiff((frac(4)+1:frac(5))', [left; i51; i52; i53; i54; i55]);

temp3 = a(i52)*a(i52+1)/(a(i52)+a(i52+1)); r(i52,i52) = r(i52,i52) + temp3; r(i52,i52+1) = r(i52,i52+1) - temp3;
temp4 = a(i53)*a(i53-1)/(a(i53)+a(i53-1)); r(i53,i53) = r(i53,i53) + temp4; r(i53,i53-1) = r(i53,i53-1) - temp4;
temp5 = a(i54)*a(i54+1)/(a(i54)+a(i54+1)); r(i54,i54) = r(i54,i54) + temp5; r(i54,i54+1) = r(i54,i54+1) - temp5;
temp6 = a(i55)*a(i55-1)/(a(i55)+a(i55-1)); r(i55,i55) = r(i55,i55) + temp6; r(i55,i55-1) = r(i55,i55-1) - temp6; 
for i = 1:length(others)
    j = others(i); temp7 = a(j)*a(j-1)/(a(j)+a(j-1)); temp8 = a(j)*a(j+1)/(a(j)+a(j+1));
    r(j,j) = r(j,j) + temp7 + temp8; r(j,j-1) = r(j,j-1) - temp7; r(j,j+1) = r(j,j+1) - temp8;
end
clear others temp3 temp4 temp5 temp6 left temp7 temp8

%% Fracture 4
left = frac(3)+1; i41 = frac(3)+2*2^refineMeshStep; i42 = i41 + 1; righ = frac(4); others = setdiff((frac(3)+1:frac(4))', [left; i41; i42; righ]);
temp1 = a(left)*a(left+1)/(a(left)+a(left+1)); r(left,left) = r(left,left) + temp1; r(left,left+1) = r(left,left+1) - temp1;
temp4 = a(righ)*a(righ-1)/(a(righ)+a(righ-1)); r(righ,righ) = r(righ,righ) + temp4; r(righ,righ-1) = r(righ,righ-1) - temp4;
temp2 = a(i41)*a(i41-1)/(a(i41)+a(i41-1)); r(i41,i41) = r(i41,i41) + temp2; r(i41,i41-1) = r(i41,i41-1) - temp2;
temp3 = a(i42)*a(i42+1)/(a(i42)+a(i42+1)); r(i42,i42) = r(i42,i42) + temp3; r(i42,i42+1) = r(i42,i42+1) - temp3;
for i = 1:length(others)
    j = others(i); temp5 = a(j)*a(j-1)/(a(j)+a(j-1)); temp6 = a(j)*a(j+1)/(a(j)+a(j+1));
    r(j,j) = r(j,j) + temp5 + temp6; r(j,j-1) = r(j,j-1) - temp5; r(j,j+1) = r(j,j+1) - temp6;
end
clear left righ others temp1 temp2 temp3 temp4 temp5 temp6

%% Fracture 10
if refineMeshStep==0
    i101 = frac(9)+1; i102 = i101+1; i103 = frac(9)+7; i104 = frac(10); left = []; righ = [];
else
    i101 = frac(9)+2^refineMeshStep; i102 = i101+1; i103 = frac(9)+7*2^refineMeshStep; i104 = i103+1;
    left = frac(9)+1; righ = frac(10);
    temp1 = a(left)*a(left+1)/(a(left)+a(left+1)); r(left,left) = r(left,left) + temp1; r(left,left+1) = r(left,left+1) - temp1;
    temp4 = a(righ)*a(righ-1)/(a(righ)+a(righ-1)); r(righ,righ) = r(righ,righ) + temp4; r(righ,righ-1) = r(righ,righ-1) - temp4;
    temp2 = a(i101)*a(i101-1)/(a(i101)+a(i101-1)); r(i101,i101) = r(i101,i101) + temp2; r(i101,i101-1) = r(i101,i101-1) - temp2;
    temp3 = a(i104)*a(i104+1)/(a(i104)+a(i104+1)); r(i104,i104) = r(i104,i104) + temp3; r(i104,i104+1) = r(i104,i104+1) - temp3;
    clear temp1 temp4 temp2 temp3
end
others = setdiff((frac(9)+1:frac(10))', [left; i101; i102; i103; i104; righ]);
temp5 = a(i102)*a(i102+1)/(a(i102)+a(i102+1)); r(i102,i102) = r(i102,i102) + temp5; r(i102,i102+1) = r(i102,i102+1) - temp5;
temp6 = a(i103)*a(i103-1)/(a(i103)+a(i103-1)); r(i103,i103) = r(i103,i103) + temp6; r(i103,i103-1) = r(i103,i103-1) - temp6;
for i = 1:length(others)
    j = others(i); temp7 = a(j)*a(j-1)/(a(j)+a(j-1)); temp8 = a(j)*a(j+1)/(a(j)+a(j+1));
    r(j,j) = r(j,j) + temp7 + temp8; r(j,j-1) = r(j,j-1) - temp7; r(j,j+1) = r(j,j+1) - temp8;
end    
clear left righ others temp5 temp6 temp7 temp8

%% Fracture 8
if refineMeshStep==0
    i81 = frac(7)+1; i82 = i81+1; left = [];
else
    i81 = frac(7)+2^refineMeshStep; i82 = i81+1; left = frac(7)+1;
    temp1 = a(left)*a(left+1)/(a(left)+a(left+1)); r(left,left) = r(left,left) + temp1; r(left,left+1) = r(left,left+1) - temp1;
    temp2 = a(i81)*a(i81-1)/(a(i81)+a(i81-1)); r(i81,i81) = r(i81,i81) + temp2; r(i81,i81-1) = r(i81,i81-1) - temp2;
    clear temp1 temp2
end
i83 = frac(7)+10*2^refineMeshStep; i84 = i83+1; righ = frac(8);
others = setdiff((frac(7)+1:frac(8))', [left; i81; i82; i83; i84; righ]);
temp6 = a(i82)*a(i82+1)/(a(i82)+a(i82+1)); r(i82,i82) = r(i82,i82) + temp6; r(i82,i82+1) = r(i82,i82+1) - temp6;
temp3 = a(i83)*a(i83-1)/(a(i83)+a(i83-1)); r(i83,i83) = r(i83,i83) + temp3; r(i83,i83-1) = r(i83,i83-1) - temp3;
temp4 = a(i84)*a(i84+1)/(a(i84)+a(i84+1)); r(i84,i84) = r(i84,i84) + temp4; r(i84,i84+1) = r(i84,i84+1) - temp4;
temp5 = a(righ)*a(righ-1)/(a(righ)+a(righ-1)); r(righ,righ) = r(righ,righ) + temp5; r(righ,righ-1) = r(righ,righ-1) - temp5;
for i = 1:length(others)
    j = others(i); temp7 = a(j)*a(j-1)/(a(j)+a(j-1)); temp8 = a(j)*a(j+1)/(a(j)+a(j+1));
    r(j,j) = r(j,j) + temp7 + temp8; r(j,j-1) = r(j,j-1) - temp7; r(j,j+1) = r(j,j+1) - temp8;
end 
clear others left righ temp3 temp4 temp5 temp6 temp7 temp8


%% Deal with intersect points
angle = [5.119253873298953e-01     1.610655443491531e+00     2.123410727191778e+00
         7.600410772429374e-01     8.564696224717313e-01     8.520317395310949e-01];
inte = [i11 i22 i12 i21;
        i101 i41 i102 i42;
        i103 i81 i104 i82;
        i51 i84 i52 i83;
        i53 i72 i54 i71];
for i = 1:size(inte,1)
    an = angle(i);
    angle2 = [0 pi-an pi an;
              pi-an 0 an pi;
              pi an 0 pi-an;
              an pi pi-an 0];
    allTran = a(inte(i,1)) + a(inte(i,2)) + a(inte(i,3)) + a(inte(i,4));
    for k = 1:4
        temp = setdiff(1:4, k);
        r(inte(i,k),inte(i,k)) = r(inte(i,k),inte(i,k)) + a(inte(i,k))*a(inte(i,temp(1)))/allTran*cos(angle2(k,temp(1))/2) + ...
                                                          a(inte(i,k))*a(inte(i,temp(2)))/allTran*cos(angle2(k,temp(2))/2) + ...
                                                          a(inte(i,k))*a(inte(i,temp(3)))/allTran*cos(angle2(k,temp(3))/2);
        r(inte(i,k),inte(i,temp(1))) = r(inte(i,k),inte(i,temp(1))) - ...
                                       a(inte(i,k))*a(inte(i,temp(1)))/allTran*cos(angle2(k,temp(1))/2);
        r(inte(i,k),inte(i,temp(2))) = r(inte(i,k),inte(i,temp(2))) - ...
                                       a(inte(i,k))*a(inte(i,temp(2)))/allTran*cos(angle2(k,temp(2))/2);
        r(inte(i,k),inte(i,temp(3))) = r(inte(i,k),inte(i,temp(3))) - ...
                                       a(inte(i,k))*a(inte(i,temp(3)))/allTran*cos(angle2(k,temp(3))/2);
    end
end  

temp1 = a(i62)*a(i55)/(a(i62)+a(i55))*cos(angle(6)/2);
r(i62,i62) = r(i62,i62) + temp1; r(i55,i55) = r(i55,i55) + temp1;
r(i62,i55) = r(i62,i55) - temp1; r(i55,i62) = r(i55,i62) - temp1;

end