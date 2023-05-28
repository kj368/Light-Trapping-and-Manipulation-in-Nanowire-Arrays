%% CODE FOR FIGURE 5 %%

% INITIALIZE MATLAB
close all;
clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a, Y] = meshgrid(0.0001:0.01:0.5, 1:350);
d = 5;

%% Calculate C1 & C2
C1 = C_parallel(a,d);
C2 = C_perpendicular(a,d);

%% Calculate PER 
A = (C1.*(2.*pi./Y).*d).^2;
B = (C2.*(2.*pi./Y).*d).^2;
T1 = A./(1+A);
T2 = 1./(1+B);
PER = 10.*log10(T2./T1);

%% PLOTS
pcolor(a*2,Y,PER);
hold on;
[C, hContour] = contour(a*2,Y,PER,0:5:40, 'LineColor','R', 'LineWidth',2);
clabel(C, hContour, 'FontSize', 25, 'Color', 'w', 'LabelSpacing', 2000);
shading interp;
set(h, 'ylim', [0,400])
ylabel('Wavelength (\mum)')
xlabel('Nanowire diameter (\mum)')
set(gca,'fontsize',25);

%% FUNCTIONS
function C1 = C_parallel(a,d)
    C1 = (1./(2.*pi)).*log(2.*(1-cos(2.*pi.*a./d)));
end
function C2 = C_perpendicular(a,d)
    C2 = (1./(2.*pi)).*(2.*(1-cosh(2.*pi.*a./d)));
end
