function vel = LSC_DGS_Stokes_typeMHD(vel, rhs, A, B, Lap_p)

n = sqrt(size(B,1)); dofu = n*(n+1); dofp = n^2;
u = vel(1:2*dofu); p = vel(2*dofu+1:end);
f = rhs(1:2*dofu); g = rhs(2*dofu+1:end);

%% Relax momentum equation
du = zeros(2*dofu,1);
du = du + tril(A) \ (f - A*u - B'*p - A*du);
du = du + triu(A) \ (f - A*u - B'*p - A*du);
u = u + du;

%% Relax transformed continuity equation
dq = zeros(dofp,1);
dq = dq + tril(Lap_p) \ (g - B*u - Lap_p*dq);
dq = dq + triu(Lap_p) \ (g - B*u - Lap_p*dq);

%% Transform the correction back to the original variables
u = u + B'*dq;
temp = zeros(dofp,1);
temp = temp + tril(Lap_p) \ (B*A*B'*dq - Lap_p*temp);
temp = temp + triu(Lap_p) \ (B*A*B'*dq - Lap_p*temp);
p = p - temp;
p = p - mean(p);

vel = [u; p];

end