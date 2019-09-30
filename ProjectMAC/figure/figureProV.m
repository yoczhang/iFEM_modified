function figureProV

n = 10; h = 2 / n; H = 2*h; N = n/2;
[vx,vy] = meshgrid(-1+h/2:h:1-h/2,1:-h:-1);
[vx2,vy2] = meshgrid(-1+H/2:H:1-H/2,1:-H:-1);
[node,elem] = squarequadmesh([-1,1,-1,1],h);
[node2,elem2] = squarequadmesh([-1,1,-1,1],H);


%% u
figure; title('Coarse grid'); showmesh(node,elem); hold on;
% ,'Edgecolor','r'

for i = 1:N+1
    for j = 1:N
        plot(vx2(i,j), vy2(i,j), 'ko', 'LineWidth', 10, 'MarkerSize', 2, 'MarkerFaceColor', 'k'); hold on;
    end
end
hold off
saveas(gcf,'106.fig');

%% Coarse grid
open 106.fig
title('1');
hold on
for i = 2:N
    for j = 1:N-1
        plot(vx(2*i-1,2*j), vy(2*i-1,2*j), 'ro', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r'); hold on;
        plot(vx(2*i-1,2*j+1), vy(2*i-1,2*j+1), 'ro', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r'); hold on;
    end
end
hold off
saveas(gcf,'107.fig');


open 107.fig
title('2');
hold on
for i = 2:N+1
    for j = 1:N-1
        plot(vx(2*i-2,2*j), vy(2*i-2,2*j), 'bo', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
        plot(vx(2*i-2,2*j+1), vy(2*i-2,2*j+1), 'bo', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
    end
end
hold off 
saveas(gcf,'108.fig');

open 108.fig
title('3');
hold on
for i = 2:N
    plot(vx(2*i-1,1), vy(2*i-1,1), 'bo', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
    plot(vx(2*i-1,end), vy(2*i-1,end), 'bo', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
end
hold off
saveas(gcf,'109.fig');

open 109.fig
title('4');
hold on
for i = 2:N+1
    plot(vx(2*i-2,1), vy(2*i-2,1), 'ro', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r'); hold on;
    plot(vx(2*i-2,end), vy(2*i-2,end), 'ro', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r'); hold on;
end
hold off
saveas(gcf,'110.fig');
end