%% CODE FOR FIGURE 3 %%

% INITIALIZE MATLAB
close all;
clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X, Y] = meshgrid(1:20, 1:350);
Data = 1./(1+((X./Y).*log(2-2.*cos(0.2.*pi./X))).^2);
pcolor(Data)
hold on 
[C, hContour] = contour(X,Y,Data,0:0.1:0.9, 'LineColor','R', 'LineWidth',2);
clabel(C, hContour, 'FontSize', 25, 'Color', 'w', 'LabelSpacing', 2000);
colorbar;
shading interp;
ylabel('Wavelength (\mum)')
xlabel('Grid period (\mum)')
set(gca,'fontsize',25);
hContour = gca;
hContourChildren = hContour.Children;
for i = 1:length(hContourChildren)
    if strcmp(hContourChildren(i).Type, 'text')
        hContourChildren(i).FontSize = 1; % Set desired font size
    end
end