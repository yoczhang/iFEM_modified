function figureIndex

n = 4; h = 2 / n;
[ux,uy] = meshgrid(-1:h:1,1-h/2:-h:-1+h/2);
[vx,vy] = meshgrid(-1+h/2:h:1-h/2,1:-h:-1);
[px,py] = meshgrid(-1+h/2:h:1-h/2,1-h/2:-h:-1+h/2);

[node,elem] = squarequadmesh([-1,1,-1,1],h);

%% p, u, v
figure; title('u(red), v(blue), p(black)'); showmesh(node,elem); hold on;
for i = 1:n
    for j = 1:n+1
        plot(ux(i,j), uy(i,j), 'ro', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r'); hold on;
    end
end
for i = 1:n+1
    for j = 1:n
        plot(vx(i,j), vy(i,j), 'bo', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b'); hold on;
    end
end
for i = 1:n
    for j = 1:n
        plot(px(i,j), py(i,j), 'ko', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'k'); hold on;
    end
end


%% p
figure; title('p(black)'); showmesh(node,elem); hold on;
for i = 1:n
    for j = 1:n
        plot(px(i,j), py(i,j), 'ko', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'k'); hold on;
    end
end

%% u
figure; title('u(red)'); showmesh(node,elem); hold on;
for i = 1:n
    for j = 1:n+1
        plot(ux(i,j), uy(i,j), 'ro', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r'); hold on;
    end
end

%% v
figure; title('v(blue)'); showmesh(node,elem); hold on;
for i = 1:n+1
    for j = 1:n
        plot(vx(i,j), vy(i,j), 'bo', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b'); hold on;
    end
end

end