%% CODE FOR FIGURE 2 %%

% INITIALIZE MATLAB
close all;
clc;
clear all;
%K=19, a =1;
[x, y] = meshgrid(-2.5:0.4:2,-2:0.4:2);
r= sqrt(x.^2 + y.^2) >1;

U = 10.*(1+ (x.^2 - y.^2)./((x.^2 + y.^2).^2));
V = 10.*((2.*x.*y)./((x.^2 + y.^2).^2));
U = U .* r;
V = V .* r;
quiver(x,y,U,V, 'r')
hold on;
p = nsidedpoly(1000, 'Center', [0, 0], 'Radius', 1);
plot(p, 'FaceColor', 'g', 'FaceAlpha',.1)
axis equal
xlabel('x');
ylabel('y');
set(gca,'fontsize',25);











