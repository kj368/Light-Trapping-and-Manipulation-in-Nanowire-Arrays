%% CODE FOR FIGURE 1 %%

% INITIALIZE MATLAB
close all;
clc;
clear all;
%K=19, a =1;
[x, y] = meshgrid(-2:0.4:2.0,-2.0:0.4:2.0);
r= sqrt(x.^2 + y.^2) <1;
U = 0.1;
V = 0;
U = U .* r;
V = V .* r;
quiver(x,y,U,V,0, 'r');
hold on
r0= sqrt(x.^2 + y.^2) >1;
U0 = (1+ (0.9).*(x.^2 - y.^2)./((x.^2 + y.^2).^2));
V0 = (((1.8).*x.*y)./((x.^2 + y.^2).^2));
U0 = U0 .* r0;
V0 = V0 .* r0;
quiver(x,y,U0,V0, 'r');
hold on;
p = nsidedpoly(1000, 'Center', [0, 0], 'Radius', 1);
plot(p, 'FaceColor', 'g', 'FaceAlpha',.1)
axis equal
xlabel('x');
ylabel('y');
set(gca,'fontsize',25);
