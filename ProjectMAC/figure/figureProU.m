function figureProU

n = 8; h = 2 / n; H = 2*h; N = n/2;
[ux,uy] = meshgrid(-1:h:1,1-h/2:-h:-1+h/2);
[ux2,uy2] = meshgrid(-1:H:1,1-H/2:-H:-1+H/2);
[node,elem] = squarequadmesh([-1,1,-1,1],h);
[node2,elem2] = squarequadmesh([-1,1,-1,1],H);

%% u
figure; title('Coarse grid'); showmesh(node,elem); hold on;
% ,'Edgecolor','r'

for i = 1:N
    for j = 1:N+1
        plot(ux2(i,j), uy2(i,j), 'ko', 'LineWidth', 10, 'MarkerSize', 2, 'MarkerFaceColor', 'k'); hold on;
    end
end
hold off
saveas(gcf,'101.fig');

%% Coarse grid 
open 101.fig
title('1');
hold on
for i = 1:N-1
    for j = 2:N
        plot(ux(2*i,2*j-1), uy(2*i,2*j-1), 'ro', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r'); hold on;
        plot(ux(2*i+1,2*j-1), uy(2*i+1,2*j-1), 'ro', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r'); hold on;
    end
end
hold off
saveas(gcf,'102.fig');

open 102.fig
title('2');
hold on
for i = 1:N-1
    for j = 2:N+1
        plot(ux(2*i,2*j-2), uy(2*i,2*j-2), 'bo', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
        plot(ux(2*i+1,2*j-2), uy(2*i+1,2*j-2), 'bo', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
    end
end
hold off
saveas(gcf,'103.fig');

open 103.fig
title('3');
hold on
for j = 2:N
    plot(ux(1,2*j-1), uy(1,2*j-1), 'bo', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
    plot(ux(end,2*j-1), uy(end,2*j-1), 'bo', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
end
hold off
saveas(gcf,'104.fig');

open 104.fig
title('4');
hold on
for j = 2:N+1
    plot(ux(1,2*j-2), uy(1,2*j-2), 'ro', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r'); hold on;
    plot(ux(end,2*j-2), uy(end,2*j-2), 'ro', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r'); hold on;
end
hold off
saveas(gcf,'105.fig');
end