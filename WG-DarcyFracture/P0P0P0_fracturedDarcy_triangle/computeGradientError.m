function value = computeGradientError(node, elem, elem2edge, area, error0, errorb)

l1 = sum((node(elem(:,3),:) - node(elem(:,2),:)).^2, 2); 
l2 = sum((node(elem(:,1),:) - node(elem(:,3),:)).^2, 2); 
l3 = sum((node(elem(:,2),:) - node(elem(:,1),:)).^2, 2); 
l12 = l1 + l2;
l23 = l2 + l3;
l13 = l3 + l1; 
l123 = l1 + l2 + l3;

v0 = error0; 
v1 = errorb(elem2edge(:,1));
v2 = errorb(elem2edge(:,2));
v3 = errorb(elem2edge(:,3));

cof1 = (32*area.^2.*(v1+v2+v3-3*v0) + 2*l1.*l123.*v1 + (l3-l12).*l123.*v2 + (l2-l13).*l123.*v3) ./ (4*area.^2.*l123);
cof2 = (32*area.^2.*(v1+v2+v3-3*v0) + 2*l2.*l123.*v2 + (l3-l12).*l123.*v1 + (l1-l23).*l123.*v3) ./ (4*area.^2.*l123);
cof3 = (32*area.^2.*(v1+v2+v3-3*v0) + 2*l3.*l123.*v3 + (l2-l13).*l123.*v1 + (l1-l23).*l123.*v2) ./ (4*area.^2.*l123);

node1 = node(elem(:,1),:); node2 = node(elem(:,2),:); node3 = node(elem(:,3),:);
mid1 = (node2 + node3)/2; mid2 = (node3 + node1)/2; mid3 = (node1 + node2)/2;

value = area/3 .* (temp(mid1, cof1, cof2, cof3, node1, node2, node3) + ...
                   temp(mid2, cof1, cof2, cof3, node1, node2, node3) + ...
                   temp(mid3, cof1, cof2, cof3, node1, node2, node3));
value = sqrt(sum(value,1));

    function r = temp(mid, co1, co2, co3, nv1, nv2, nv3)
        x = mid(:,1); y = mid(:,2);
        r = (co1.*(x-nv1(:,1))+co2.*(x-nv2(:,1))+co3.*(x-nv3(:,1))).^2 + ...
            (co1.*(y-nv1(:,2))+co2.*(y-nv2(:,2))+co3.*(y-nv3(:,2))).^2;
    end
end