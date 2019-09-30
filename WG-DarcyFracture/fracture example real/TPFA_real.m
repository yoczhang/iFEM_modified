function [r, b] = TPFA_real(NGamma, lGamma, Kf, fracLength, line, node, leftFractureEdge, rightFractureEdge)

co = lGamma*Kf; alpha = 2*co./fracLength;
r = sparse(NGamma, NGamma); b = zeros(NGamma,1);

for i = 1:NGamma
    bp = line(i,1); ep = line(i,2);
    [flagBp, unused] = find(line==bp);
    [flagEp, unused] = find(line==ep);
    flagBp = setdiff(flagBp, i); flagEp = setdiff(flagEp, i);
    
    % length(flagBp==0) denotes no intersection edges
    if length(flagBp)==1 % denotes only one intersection edge
        j = flagBp;
        r(i,i) = r(i,i) + alpha(i)*alpha(j) / (alpha(i)+alpha(j));
        r(i,j) = r(i,j) - alpha(i)*alpha(j) / (alpha(i)+alpha(j));
    elseif length(flagBp)==3 % denotes three intersection edges
        j1 = flagBp(1); j2 = flagBp(2); j3 = flagBp(3);
        tranSum = alpha(i) + alpha(j1) + alpha(j2) + alpha(j3);
        vi  = node(line(i,1),:)  - node(line(i,2),:);
        vj1 = node(line(j1,1),:) - node(line(j1,2),:);
        vj2 = node(line(j2,1),:) - node(line(j2,2),:);
        vj3 = node(line(j3,1),:) - node(line(j3,2),:);
        
        angle1 = acos(vi*vj1'/(norm(vi)*norm(vj1)));
        angle2 = acos(vi*vj2'/(norm(vi)*norm(vj2)));
        angle3 = acos(vi*vj3'/(norm(vi)*norm(vj3)));
        
        r(i,i)  = r(i,i)  + alpha(i)*alpha(j1) / tranSum * cos(angle1/2) + ...
                            alpha(i)*alpha(j2) / tranSum * cos(angle2/2) + ...
                            alpha(i)*alpha(j3) / tranSum * cos(angle3/2);
        r(i,j1) = r(i,j1) - alpha(i)*alpha(j1) / tranSum * cos(angle1/2);
        r(i,j2) = r(i,j2) - alpha(i)*alpha(j2) / tranSum * cos(angle2/2);
        r(i,j3) = r(i,j3) - alpha(i)*alpha(j3) / tranSum * cos(angle3/2);
    end
    
    if length(flagEp)==1
        j = flagEp;
        r(i,i) = r(i,i) + alpha(i)*alpha(j) / (alpha(i)+alpha(j));
        r(i,j) = r(i,j) - alpha(i)*alpha(j) / (alpha(i)+alpha(j));
    elseif length(flagEp)==3
        j1 = flagEp(1); j2 = flagEp(2); j3 = flagEp(3);
        tranSum = alpha(i) + alpha(j1) + alpha(j2) + alpha(j3);
        vi  = node(line(i,1),:)  - node(line(i,2),:);
        vj1 = node(line(j1,1),:) - node(line(j1,2),:);
        vj2 = node(line(j2,1),:) - node(line(j2,2),:);
        vj3 = node(line(j3,1),:) - node(line(j3,2),:);
        
        angle1 = acos(vi*vj1'/(norm(vi)*norm(vj1)));
        angle2 = acos(vi*vj2'/(norm(vi)*norm(vj2)));
        angle3 = acos(vi*vj3'/(norm(vi)*norm(vj3)));
        
        r(i,i)  = r(i,i)  + alpha(i)*alpha(j1) / tranSum * cos(angle1/2) + ...
                            alpha(i)*alpha(j2) / tranSum * cos(angle2/2) + ...
                            alpha(i)*alpha(j3) / tranSum * cos(angle3/2);
        r(i,j1) = r(i,j1) - alpha(i)*alpha(j1) / tranSum * cos(angle1/2);
        r(i,j2) = r(i,j2) - alpha(i)*alpha(j2) / tranSum * cos(angle2/2);
        r(i,j3) = r(i,j3) - alpha(i)*alpha(j3) / tranSum * cos(angle3/2);
    end
end

% left and right sides, all have a boundary points,
% since gD = 0 on right side, we don't deal with it,
% on left side, gD = 101325, we must find the number
% of edge on the left side in all fractured edges 
r(leftFractureEdge,leftFractureEdge)   = r(leftFractureEdge,leftFractureEdge)   + alpha(leftFractureEdge);
r(rightFractureEdge,rightFractureEdge) = r(rightFractureEdge,rightFractureEdge) + alpha(rightFractureEdge);
b(leftFractureEdge) = alpha(leftFractureEdge)*1013250;
end