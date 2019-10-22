function check_results(uh, vh, uI, vI)


n = min(size(uh));
uh = uh(:); vh = vh(:); 
ue = uI - uh; ve = vI - vh;

%- check uh
nx = n+1;  ny = n;
LEF = 1:ny;   RIG = ny*(nx-1)+1:ny*nx;
TOP = ny+1:ny:ny*(nx-2)+1; BOT = ny*2:ny:ny*(nx-1);
a = false;  
a(LEF) = true; 
a(RIG) = true; 
inteNodes = find(~a);

ue(LEF)
ue(RIG)

%- check vh
nx = n;  ny = n+1;
TOP = 1:ny:ny*(nx-1)+1;  BOT = ny:ny:ny*nx;
LEF = 2:ny-1;  RIG = ny*(nx-1)+2:ny*nx-1;
a = false;  
a(TOP) = true;  
a(BOT) = true;  
inteNodes = find(~a);

ve(TOP)
ve(BOT)

end