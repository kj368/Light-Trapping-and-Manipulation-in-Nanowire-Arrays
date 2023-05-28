%% CODE FOR FIGURE 6 %%

% INITIALIZE MATLAB
close all;
clc;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The dilectric is considered to only have a permittivity and has a
% permeabilty of 1. However the magnetic Field is still included in the 
% simultion for completeness but can be ignored.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NGP = Number of grid points
% GR = Grid resolution 
% TS = Time step
% NI = Number of iterations
% er1 = Electric permittivity Of dilectric
% er2 = Electric permittivity everywhere else
% epsilon_array =  Material electric permittivity array using er1
% mu_array = Material magnetic permeability array
% T = Thickness of the dilectric
% epsilon_max = Maximum Electric permittivity
% mu_r = Local magnetic permeability at time instance
% EF = Electric field
% HF = Magnetic field
% EFUC = Electric field update coefficient
% HFUC = Magnetic field update coefficient
% MAX_freq = Maximum frequency 
% MAX_n = maximum refractive index
% R_lambda = Resolution due to wavelength
% NSPC = Number of cells before and after device
% MIN_dim = Absolute minimum dimension 
% MIN_lambda = Minimum wavelength 
% GR_con1 = Condition one for grid resolution
% GR_con2 = Condition two for grid resolution
% NUM_diel = Number of cells for dielectric
% IND_start = Start cell index for dilectric
% IND_end = End cell index for dilectric
% Z = Z-axis
% n_bound = Refractive index at the boundary
% POS_source = Position of the guassian source
% tau = Duration of the guassian source
% t0 = Offset of the guassian source
% t = Time array
% EF_source = Electric field guassian source
% HF_source = Magnetic field guassian source
% E1 = Electric field at first time step
% E2 = Electric field at second time step
% H1 = Magnetic field at first time step
% H2 = Magnetic field at second time step
% A = Amplitude of Magetic field (slightly varied from electric field)
% s = Delay of Magetic field (slightly varied from electric field)
% n_source = Refractive index of postion of the source injection
% NUM_freq = Number of frequency points oin fourier transform
% FREQ_array = Frequency array for fourier transform
% K = Kernal for each frequency 
% FT_ref = Array of Reflection fourier transform
% FT_trn = Array of Transmission fourier transform
% FT_source = Array of Source fourier transform
% REF = Reflectivity
% TRN = Transmission
% CON = Control 




% CONSTANTS
c0 = 3*10^8;
e0 = 8.85*10^-12;
u0 = 1.26*10^-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SOURCE PARAMETERS 
fmax = 5*(10^9);
NFREQ = 1000;
FREQ = linspace(0, fmax, NFREQ);

% DEVICE PROPERTIES
er1 = 12.0;
er2 = 1.0;
L = 3.0*10^-2;

% GRID PARAMETERS
ermax = max(er1, er2);
nmax = sqrt(ermax);
NRES_LAM = 20;
NRES_D = 2;
NSPC = [100, 100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calulate grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial resolution
lam_min = c0/fmax/nmax;
dz1 = lam_min/NRES_LAM;
dz2 = L/NRES_D;
dz = min([dz2, dz1]);

% CRITICAL DIMENSIONS
nz = ceil(L/dz);
dz = L/nz;

% NUMBER OF GRID CELL
Nz= nz + sum(NSPC) + 3;  % Buffer 3 introduced


% FIELD POSITION
za = [0:Nz-1]*dz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEVICE PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITILIZE RELATIVE MATERIAL PROPERTIES
ERxx = er2*ones(1,Nz);          % Relative permittivity
URyy = ones(1,Nz);              % Relative permeability

% POSITION OF SLAB
nz1 = 2 + NSPC(1) + 1;
nz2 = nz1 + round(L/dz) - 1;

% PROPERTIES OF DIELECTRIC
ERxx(nz1:nz2) = er1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TIME STEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% REFRACTIVE INDEX AT BOUNDARY
nbc = sqrt(URyy(1)*ERxx(1));

% Timestep for PABC
dt = nbc*dz/(2*c0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GAUSSIAN SOURCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_src = 20;
tau = 0.5/fmax;
t0 = 3*tau;
N = 3000;                           % Number of time steps
t = [0:N-1]*dt;
nsrc = sqrt(ERxx(k_src)*URyy(k_src));
s = 0.5*nsrc*dz/c0 -dt/2;           % delay for Hysrc
A = sqrt(ERxx(k_src)/URyy(k_src));  % amplitude of Hysrc
ExSrc = exp(-((t-t0)/tau).^2);
HySrc = A*exp(-((t-t0+s)/tau).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FOURIER TRANSFORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kernal
K = exp((-1i*2*pi*dt)*FREQ);
% Initialise FT
ExR = zeros(1,NFREQ);
ExT = zeros(1,NFREQ);
SRC = zeros(1,NFREQ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIELD INITIALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ex = zeros(1,Nz);
Hy = zeros(1,Nz);

% INITIALISE BOUNDARY
E1 = 0;  E2 = 0;
H1 = 0;  H2 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UPDATE COEFFICEINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mHy = -(c0*dt/dz)./URyy;        
mEx = -(c0*dt/dz)./ERxx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:N
    % Ex at boundary
    E2 = E1;
    E1 = Ex(Nz);

    % UPDATE Ex
    Ex(1) = Ex(1) + mEx(1)*(Hy(1) - H2);
    for k = 2:Nz
        Ex(k) = Ex(k) + mEx(k)*(Hy(k)- Hy(k-1));
    end

    % ADD ELECTRIC SOURCE & TF/SF correction
    Ex(k_src) = Ex(k_src) -mEx(k_src)*HySrc(n);
    
    % Hy at boundary
    H2 = H1;
    H1 = Hy(1);
  
    % UPDATE Hy
    for k = 1:Nz-1 
        Hy(k) = Hy(k) + mHy(k)*(Ex(k+1)- Ex(k));
    end
    Hy(Nz) = Hy(Nz) +mHy(Nz)*(E2 - Ex(Nz));
    
    % ADD MAGNETIC SOURCE & TF/SF correction
    Hy(k_src-1) = Hy(k_src-1) - mHy(k_src-1)*ExSrc(n);

    % Update FT 
    for nf = 1:NFREQ
        ExR(nf) = ExR(nf) + (K(nf)^n)*Ex(1);
        ExT(nf) = ExT(nf) + (K(nf)^n)*Ex(Nz);
        SRC(nf) = SRC(nf) + (K(nf)^n)*ExSrc(n);
    end

    % Visualise
    if ~mod(n,10)
        % Calculate Spectra
        REF = abs(ExR./SRC).^2;
        TRN = abs(ExT./SRC).^2;
        CON = REF+TRN;

        % Show field
        clf;
        subplot(212);
        plot(za,Ex,'-b');
        hold on;
        plot(za,Hy,'-r');
        hold on;
        x = [0.0877292, 0.11668838, 0.11668838, 0.0877292,];
        y = [-2, -2, 2, 2];
        patch(x,y,'g','FaceAlpha',.2,'EdgeAlpha',.4)
        hold on;
        x = [0.0170348, 0.0170348];
        y = [-2, 2];
        patch(x,y,'r')
        xlim([za(1), za(Nz)]);
        ylim([-1.8,1.8]);
        ylabel('amplitude')
        xlabel('z-axis')
        title(['step ', num2str(n), ' of ', num2str(N)]);
        set(gca,'fontsize',32);
        
        % Show Spectra
        subplot(211);
        plot(FREQ,REF,'-r', 'LineWidth',3);
        hold on;
        plot(FREQ,TRN,'-b', 'LineWidth',3);
        plot(FREQ,CON,':k', 'LineWidth',3);
        hold off;
        xlim([0, fmax]);
        ylim([0,1.3]);
        ylabel('ratio')
        xlabel('Frequency')
        title('spectra')
        set(gca,'fontsize',32);
        drawnow;
    end
end

