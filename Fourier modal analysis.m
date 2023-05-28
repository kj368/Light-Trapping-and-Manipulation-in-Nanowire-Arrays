% INITILISE MATLAB
close all;
clc;
clear all;

% Define the wavelength and the wavevector of the incident plane wave
lambda = 3e-4;
k0 = 2*pi/lambda;

% Define the dielectric constant of the wire
epsilon_wire = 4e12;

% Define the radius of the wire
a = 0.1*lambda;

% Define the scattering matrix
S = zeros(2,2);

% Define the x-component of the incident electric field
Einc = [1; 0];

% Define the x-component of the scattered electric field
Esc = [0; 0];

% Calculate the scattering matrix
S(1,1) = (k0*a*epsilon_wire*besselj(1,k0*a))/(k0*a*epsilon_wire*besselj(1,k0*a)+besselh(1,2,k0*a));
S(1,2) = (2*besselj(1,k0*a))/(k0*a*epsilon_wire*besselj(1,k0*a)+besselh(1,2,k0*a));
S(2,1) = (2*besselj(1,k0*a))/(k0*a*epsilon_wire*besselj(1,k0*a)+besselh(1,2,k0*a));
S(2,2) = (k0*a*epsilon_wire*besselj(1,k0*a)-besselh(1,2,k0*a))/(k0*a*epsilon_wire*besselj(1,k0*a)+besselh(1,2,k0*a));

% Calculate the scattered electric field
Esc = S*Einc;

% Magnitude of the electric field outside far from the dielectric
E_out0= abs(Einc(1));

% Magnitude of the electric field outside the dielectric
E_out= abs(Esc(1));

% Calculate the magnitude of the electric field inside the dielectric
E_in = abs(Esc(2));

[x, y] = meshgrid(-3:0.01:3,-3:0.01:3);
r = sqrt(x.^2 + y.^2) < 1;
I = E_in.* r;
r = (1 < sqrt(x.^2 + y.^2)) & (sqrt(x.^2 + y.^2) < 2);
O = E_out.* r;
r = sqrt(x.^2 + y.^2) > 2;
F = E_out0.* r;
M = F + O + I;
sigma = 20; % adjust sigma value as desired
M_smooth = imgaussfilt(M, sigma);
pcolor(M_smooth);
axis tight;
colorbar;
colormap(jet)
shading flat;
set(gca,'fontsize',25);


