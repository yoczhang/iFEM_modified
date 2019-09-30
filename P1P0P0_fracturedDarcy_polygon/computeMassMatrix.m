function M = computeMassMatrix(localEdgeLength, left, right, localArea, cen, hE, Normal)

% t = [-1, 0, 1]; weight = [1/3, 4/3, 1/3];
t = [0, sqrt(21)/7, -sqrt(21)/7, 1, -1];  weight = [32/45, 49/90, 49/90, 1/10, 1/10];
in = zeros(length(localEdgeLength), 5);
for k = 1:length(t)
    pxy = (right - left)/2 * t(k) + (right + left)/2;
    in = in + weight(k) * repmat(localEdgeLength/2,1,5) .* ri(pxy);
end
in = sum(in); i1 = in(1); i2 = in(2); i3 = in(3); i4 = in(4); i5 = in(5);

M = [localArea  i1  i2;
    i1          i3  i4;
    i2          i4  i5];

    function r = ri(p)
        x = p(:,1); y = p(:,2); xE = cen(1); yE = cen(2);
        r = [1/(2*hE)*(x-xE).^2.*Normal(:,1), 1/(2*hE)*(y-yE).^2.*Normal(:,2), ...
             1/(3*hE^2)*(x-xE).^3.*Normal(:,1), 1/(2*hE^2)*(x-xE).^2.*(y-yE).*Normal(:,1), ...
             1/(3*hE^2)*(y-yE).^3.*Normal(:,2)];
    end

end