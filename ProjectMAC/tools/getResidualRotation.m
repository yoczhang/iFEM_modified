function [rh1, rh2, rh3] = getResidualRotation(uh, vh, ph, f1h, f2h, gh, A1, A2, B, U, V) 

n = size(uh,1); dof = n*(n+1);

vel = [uh(:); vh(:)]; ph = ph(:); fh = [f1h(:); f2h(:)]; gh = gh(:);
A  = [A1 U; V A2];
r = fh - A*vel - B'*ph;
rh3 = gh - B*vel;

rh1 = r(1:dof,1); rh2 = r(dof+1:2*dof,1);
rh1 = reshape(rh1, n, n+1); rh2 = reshape(rh2, n+1, n);
rh3 = reshape(rh3, n, n);

end