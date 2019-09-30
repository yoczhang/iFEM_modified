function figureIndex2(n)

h = 2 / n;
[ux,uy] = meshgrid(-1:h:1,1-h/2:-h:-1+h/2);
[vx,vy] = meshgrid(-1+h/2:h:1-h/2,1:-h:-1);
[px,py] = meshgrid(-1+h/2:h:1-h/2,1-h/2:-h:-1+h/2);
[node,elem] = squarequadmesh([-1,1,-1,1],h);

%% u, v
figure; title('u(red), v(black)'); showmesh(node,elem); hold on;
for i = 1:n
    for j = 1:n+1
        k = i + (j-1)*n;
        text(ux(i,j)+0.05, uy(i,j),int2str(k),'FontSize',12,'Color','r');
        plot(ux(i,j), uy(i,j), 'ro', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r'); hold on;
    end
end
for i = 1:n+1
    for j = 1:n
        k = i + (j-1)*(n+1);
        text(vx(i,j), vy(i,j)-0.07,int2str(k),'FontSize',12,'Color','k');
        plot(vx(i,j), vy(i,j), 'ko', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'k'); hold on;
    end
end
hold off

%% p
figure; title('P(black)'); showmesh(node,elem); hold on;
for i = 1:n
    for j = 1:n
        k = i + (j-1)*n;
        text(px(i,j)+0.05, py(i,j),int2str(k),'FontSize',12,'Color','k');
        plot(px(i,j), py(i,j), 'ko', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'k'); hold on;
    end
end
hold off
end