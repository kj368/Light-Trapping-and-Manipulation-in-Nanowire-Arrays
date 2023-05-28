% FTDT_2D_TM.m

% INITILISE MATLAB
close all;
clc;
clear all;

% CONSTANTS
e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
N0 = sqrt(u0/e0);            % Free space impedence
c0 = 2.99792458e8;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% SOURCE PARAMETERS
NFREQ = 100;                           % Number of freuency steps
freq1 = 8.5e11;                           % Initial frequency 
freq2 = 6e12;                          % Final frequency 
FREQ = linspace(freq1, freq2, NFREQ);  % Frequency array

% LOSSY SLAB PARAMETERS
f0 = mean(FREQ);                % Using the average frequency to maintain phase
lam0 = c0/f0;
d = 0.25*lam0;
er = 12;
sig = 100;

% BERENGERS PML PARAMETERS 
pml_ky = 1;
pml_ay = 1e-10;
pml_Npml = 3;
pml_R0 = 1e-8;


% GRID PARAMETERS
NRES = 20;                  % Number of points per wavelenth (the size of grid resultion)
Sx = 0.2*lam0;                % Physical width of grid
SPACER = 0.5*lam0*[1, 1];
NPML = [20, 20];            % Size of PML
nmax = sqrt(er);            % Maximum refractive index on grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE GRID RESOLUTION 
lam0_min = c0/max(FREQ);        % Shortest freesapce wavelength
dx = lam0_min/nmax/NRES;        % Resolution of grid in x
dy = lam0_min/nmax/NRES;        % Resolution of grid in y

% SNAP GRID TO CRITICAL DIMENSION
ny = ceil(d/dy);
dy = d/ny;

% COMPUTE GRID SIZE
Nx = ceil(Sx/dx);                               % Size of grid in x
Sx = Nx*dx;                                     % Physical width 
Sy = SPACER(1) + d + SPACER(2);                 % Recalculate physical width incase of rounding
Ny = NPML(1) + ceil(Sy/dy) + NPML(2);           % Size of grid in y
Sy = Ny*dy;                  

% GRID AXES
xa = [0;Nx-1]*dx;
ya = [0;Ny-1]*dy;

% 2X GRID PARAMATERS
Nx2 = 2*Nx;                 % Size of 2X grid in x
dx2 = dx/2;                 % Resolution of 2X grid in x
Ny2 = 2*Ny;                 % Size of 2X grid in y
dy2 = dy/2;                 % Resolution of 2X grid in y
xa2 = [0:Nx2-1]*dx2;
xa2 = xa2 - mean(xa2);
ya2 = [0:Ny2-1]*dy2;
[Y2, X2] = meshgrid(ya2, xa2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALISE TO VACUUM
ERzz = ones(Nx, Ny);              % Material relative permitivity 
URxx = ones(Nx, Ny);              % Material relative permeabilty 
URyy = ones(Nx, Ny);              % Material relative permeabilty
SIGzz = zeros(Nx, Ny);            % Initialise Conductivity

% BUILD SLAB
%ny1 = NPML(1) + round(SPACER(1)/dy) + 1;    % Start index
%ny2 =  ny1 + round(d/dy) - 1;                                     % Final indez
%ERzz(:, ny1:ny2) = er;
%SIGzz(:, ny1:ny2) = sig;

% SHOW DEVICE 
%subplot(1,5,1);
%imagesc(xa,ya, ERzz.');
%axis equal tight;
%colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DILECTRIC PROPERTIES
er_air = 1;
r = 3.5e-6;                             % Radius 
NP = 3;

% INITIALISE TO VACUUM
UR2 = ones(Nx2, Ny2);              % Material relative permeabilty 

% BUILD DIELECTRIC
ER2 = zeros(Nx2, Ny2);
SIG2 = zeros(Nx2, Ny2);
y1 = NPML(1)*dy + SPACER(1);
for np = 1 : NP
    y0 = y1;
    A = (X2).^2 + (Y2-y0).^2;
    ER2 = ER2 | (A < r^2);
end 

% SCALE TO REAL WORLD
ER3 = ER2;
ER2 = er_air + (er - er_air)*ER2;
SI2 = sig*ER3;

% EXTRACT 1X GRID ARRAYS
ERzz = ER2(1:2:Nx2, 1:2:Ny2);
URxx = UR2(1:2:Nx2, 2:2:Ny2);
URyy = UR2(2:2:Nx2, 1:2:Ny2);
SIGzz = SI2(1:2:Nx2, 1:2:Ny2);

% SHOW DEVICE
subplot(1, 5, 1);
imagesc(xa2, ya2, ER2.');
set(gca,'fontsize',20) 
axis equal tight off
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE SOURCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE STABLE TIME-STEP
dmin = min([dx, dy]);             % Minimum resultion of the simulation  
dt = dmin/(2*c0);                 % Time step

% CACULATE GAUSSIAN PULSE POSITION AND RECORD PLANES
ny_ref = NPML(1) + 1;             % Position where reflectance is recorded
ny_src = ny_ref + 1;              % Postion of the source along the Ny 
ny_trn = Ny - NPML(2);            % Position where transmittance is recorded

% MATERIAL PROPERTIES
ursrc = URxx(1, ny_src);
urtrn = URxx(1, ny_trn);
ersrc = ERzz(1, ny_src);
ertrn = ERzz(1, ny_trn);
nref = sqrt(ursrc*ersrc);       % Refrative index of reflection side
ntrn = sqrt(urtrn*ertrn);       % Refrative index of transmission side

% GUASSIAN PULSE PARAMETERS
tau = 0.5/max(FREQ);                   % Pulse duration
t0 = 3*tau;              % Delay of the pulse

% CACULATE NUMBER OF TIME STEPS
tprop = nmax*Sy/c0;               % Time for the wave to propagate along whole grid
t = 2*t0 + 15*tprop;               % The total time of the simulation
STEPS = ceil(t/dt);               % Number of iterations

% COMPUTE GAUSSIAN PULSE
ETAsrc = sqrt(URxx(1, ny_src)/ERzz(1, ny_src));   % Relative impedance of the source
nsrc = sqrt(URxx(1, ny_src)*ERzz(1, ny_src));     % Relative refrative index of the source
delt = 0.5*dt + 0.5*nsrc*dy/c0;                   % Times between the magnetic and electric field
t = [0:STEPS-1]*dt;
Ezsrc = exp(-((t-t0)/tau).^2);
Hxsrc = (1/ETAsrc) * exp(-((t-t0+ delt)/tau).^2);

% COMPENSATE FOR DISPERSION 
k0 = 2*pi*f0/c0;                % Free sapce wave number
f = c0*dt/sin(pi*f0*dt)*sin(k0*dy/2)/dy;            
ERzz = f*ERzz;
URxx = f*URxx;
URyy = f*URyy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE PML CONDUCTIVITY ON 2X GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITILISE PML CONDUCTIVITIES
sigy2 = zeros(Nx2, Ny2);


% ADD YLO CONDUCTIVITY
sigmax = -(pml_Npml + 1)*log(pml_R0)/(4*N0*NPML(1)*dy2);
for ny = 1 : 2*NPML(1)
    ny1 = 2*NPML(1) - ny + 1;
    sigy2(:, ny1) = sigmax*(ny/2/NPML(1))^pml_Npml;
end
% ADD YHI CONDUCTIVITY
sigmax = -(pml_Npml + 1)*log(pml_R0)/(4*N0*NPML(2)*dy2);
for ny = 1 : 2*NPML(2)
    ny1 = Ny2 - 2*NPML(2) + ny;
    sigy2(:, ny1) = sigmax*(ny/2/NPML(2))^pml_Npml;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE UPDATE COEFFICIENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M
M = c0*dt;

% CALCULATE UPDATE COEFFICIENTS FOR Bx
sigy = sigy2(1: 2: Nx2, 2: 2: Ny2);
bBxy = exp(-(sigy/pml_ky + pml_ay)*dt/e0);
cBxy = sigy./(sigy*pml_ky + pml_ay*pml_ky^2).*(bBxy - 1);


% CALCULATE UPDATE COEFFICIENTS FOR Dz
sigy = sigy2(1: 2: Nx2, 1: 2: Ny2);
bDzy =  exp(-(sigy/pml_ky + pml_ay)*dt/e0);
cDzy = sigy./(sigy*pml_ky + pml_ay*pml_ky^2).*(bDzy - 1);

% CALCULATE UPDATE COEFFICIENTS FOR Ez
Az = (dt*c0)*SIGzz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALISE FDTD TERMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALISE FILEDS TO ZERO 
Bx = zeros(Nx, Ny);
By = zeros(Nx, Ny);
Hx = zeros(Nx, Ny);
Hy = zeros(Nx, Ny);
Dz = zeros(Nx, Ny);
Ez = zeros(Nx, Ny);

% INITIALISE DERIVATIVE ARRAYS
dEzx = zeros(Nx, Ny);
dEzy = zeros(Nx, Ny);
dHxy = zeros(Nx, Ny);
dHyx = zeros(Nx, Ny);

% INITIALISE CONVOLUTIONS
psiBx_ylo = zeros(Nx, NPML(1));
psiBx_yhi = zeros(Nx, NPML(2));
psiDz_ylo = zeros(Nx, NPML(1));
psiDz_yhi = zeros(Nx, NPML(2));

% INITIALISE INTEGRATION
IEz = zeros(Nx,Ny);

% INITIALISE FOURIER TRANSFORM
K = exp(-1i*2*pi*FREQ*dt);                      % Calculate the kernals for each position on grid
EREF = zeros(Nx, NFREQ);                        % Reflected steady state fields FT
ETRN = zeros(Nx, NFREQ);                        % Transmitted steady state fields FT
SRC = zeros(1, NFREQ);                          % FT of the source

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN FDTD LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% MAIN LOOP -- ITERATE OVER TIME
%
for T = 1 : STEPS

    % CALCULATE DERIVATIVES OF E
        % dEzy
        for nx = 1 : Nx
            for ny = 1 : Ny - 1 
                dEzy(nx, ny) = (Ez(nx, ny+1) - Ez(nx, ny))/dy;
            end
            dEzy(nx, Ny) = (Ez(nx, 1) - Ez(nx, Ny))/dy;
        end
        % dEzx
        for ny = 1 : Ny
            for nx = 1 : Nx - 1 
                dEzx(nx, ny) = (Ez(nx+1, ny) - Ez(nx, ny))/dx;
            end
            dEzx(Nx, ny) = (Ez(1, ny) - Ez(Nx, ny))/dx;
        end
    % INCOPORATE TF/SF CORRECTION
    dEzy(:, ny_src-1) = dEzy(:, ny_src-1) - Ezsrc(T)/dy;
    
    % UPDATE CONVOLUTIONS FOR B FIELD UPDATES
    psiBx_ylo = bBxy(:, 1: NPML(1)).*psiBx_ylo ...
        + cBxy(:, 1: NPML(1)).*dEzy(:, 1: NPML(1));
    psiBx_yhi = bBxy(:, Ny - NPML(2) + 1: Ny).*psiBx_yhi ...
        + cBxy(:, Ny - NPML(2) + 1: Ny).*dEzy(:, Ny - NPML(2) + 1: Ny);
    
    % UPDATE B FROM E
    Bx = Bx - M*(dEzy/pml_ky);
    Bx(:, 1: NPML(1)) = Bx(:, 1: NPML(1)) ...
        - M.* psiBx_ylo;
    Bx(:, Ny - NPML(2) + 1: Ny) = Bx(:, Ny - NPML(2) + 1: Ny) ...
        - M.* psiBx_yhi;
    By = By - M*(-dEzx);

    % UPDATE H from B
    Hx = Bx./URxx;
    Hy = By./URyy;

    % CALCULATE DERIVATIVES OF H
        % dHyx
        for ny = 1 : Ny
            dHyx(1, ny) = (Hy(1, ny) - Hy(Nx, ny))/dx;
            for nx = 2 : Nx
                dHyx(nx, ny) = (Hy(nx, ny) - Hy(nx-1, ny))/dx;
            end
        end
        % dHxy
        for nx = 1 : Nx
            dHxy(nx, 1) = (Hx(nx, 1) - Hx(nx, Ny))/dy;
            for ny = 2 : Ny 
                dHxy(nx, ny) = (Hx(nx, ny) - Hx(nx, ny-1))/dy;
            end
        end
    
    % INCOPORATE TF/SF CORRECTION
    dHxy(:, ny_src) = dHxy(:, ny_src) - Hxsrc(T)/dy;


    % UPDATE CONVOLUTIONS FOR D FIELD UPDATES
    psiDz_ylo = bDzy(:, 1: NPML(1)).*psiDz_ylo ...
       + cDzy(:, 1: NPML(1)).*dHxy(:, 1: NPML(1));
    psiDz_yhi = bDzy(:, Ny - NPML(2) + 1: Ny).*psiDz_yhi ...
       + cDzy(:, Ny - NPML(2) + 1: Ny).*dHxy(:, Ny - NPML(2) + 1: Ny);
    
    % UPDATE D FROM H
    Dz = Dz + M*(dHyx - dHxy/pml_ky);
    Dz(:, 1: NPML(1)) = Dz(:, 1: NPML(1)) - M*psiDz_ylo;
    Dz(:, Ny - NPML(2) + 1: Ny) = Dz(:, Ny - NPML(2) + 1: Ny)...
        - M*psiDz_yhi;

    % UPDATE E from D
    IEz = IEz + Ez;
    Ez = (Dz - Az.*IEz)./ERzz;

    % UPDATE FOURIER TRANSFORM
    for nfreq = 1 : NFREQ
        EREF(:,nfreq) = EREF(:,nfreq) + (K(nfreq)^T)*Ez(:, ny_ref);
        ETRN(:,nfreq) = ETRN(:,nfreq) + (K(nfreq)^T)*Ez(:, ny_trn);
        SRC(nfreq) = SRC(nfreq) + (K(nfreq)^T)*Ezsrc(T);
    end
    
    % SHOW FIELDS
    d = 1;                     % decrease number to increase the darkness of the wave fronts
    if mod(T,10) == 0

        % CALCULATE REF & TRN ARRAYS
        REF = zeros(1, NFREQ);
        TRN = zeros(1, NFREQ );
        for nfreq = 1 : NFREQ

            % GET NEXT FREQUENCY 
            lam0 = c0/FREQ(nfreq);
            k0 = 2*pi/lam0;

            % CALCULATE WAVE VECTOR COMPONENTS
            m = [-floor(Nx/2) : +  floor((Nx-1)/2)].';
            kxi = -m*2*pi/Sx;
            kyinc = k0*nref;
            kyref = sqrt((k0*nref)^2 - kxi.^2);
            kytrn = sqrt((k0*ntrn)^2 - kxi.^2);

            % REFLECTION 
            eref = EREF(:,nfreq)/SRC(nfreq);
            aref = fftshift(fft(eref))/Nx;
            RDE = abs(aref).^2.*real(kyref/kyinc);
            TRN(nfreq) = sum(RDE); 

             % TRANSMISSION 
            etrn = ETRN(:,nfreq)/SRC(nfreq);
            atrn = fftshift(fft(etrn))/Nx;
            TDE = abs(atrn).^2.*real(ursrc/urtrn*kytrn/kyinc);
            REF(nfreq) = sum(TDE);
        end
        if TRN == 0
            ABSOR =  0.5* REF;
            CON = REF + ABSOR + (REF + TRN)/T;

        else
             ABSOR = 1 - (REF + TRN);
             CON = REF + TRN + ABSOR + (REF + TRN)/T;
        end
       
        % SHOW FIELD
        subplot(1, 5, 2)
        imagesc(xa, ya, Ez.');
        set(gca,'fontsize',20) 
        axis equal tight;
        colorbar;
        caxis(d*[-1, 1]);  

        % SHOW FREQUENCY RESPONSE 
        subplot(1, 5, 3:5)
        plot(FREQ*10^(-12), CON, ':k', 'LineWidth',3);
        hold on;
        plot(FREQ*10^(-12), REF, '-r','LineWidth',3);
        plot(FREQ*10^(-12), TRN, '-b','LineWidth',3);
        plot(FREQ*10^(-12), ABSOR, '-G','LineWidth',3);
        hold off;
        axis tight;
        xlim([FREQ(1), FREQ(NFREQ)]*10^(-12))
        ylim([-0.05, 1.1]);
        xlabel('FREQUENCY (THz)');
        title(['iteration ' num2str(T) ' of ' num2str(STEPS)]);
        set(gca,'fontsize',20) 
        drawnow;
    end
end
