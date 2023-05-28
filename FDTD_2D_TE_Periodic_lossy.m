% FTDT_2D_TE.m

% INITILISE MATLAB
close all;
clc;
clear all;
% CONSTANTS
e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
N0 = sqrt(u0/e0);            % Free space impedence
c0 = 2.99792458e8;

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
ER2 = ones(Nx2, Ny2);              % Material relative permitivity 
UR2 = ones(Nx2, Ny2);              % Material relative permeabilty 
SIG2 = zeros(Nx2, Ny2);            % Initialise Conductivity

% BUILD SLAB
%ny1 = 2*NPML(1) + round(SPACER(1)/dy2) + 1;    % Start index
%ny2 =  ny1 + round(d/dy2) - 1;                                     % Final indez
%ER2(:, ny1:ny2) = er;
%SIG2(:, ny1:ny2) = sig;

% DILECTRIC PROPERTIES
er_air = 1;
r = 4e-6;                             % Radius 
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
SIG2 = sig*ER3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% EXTRACT 1X GRID ARRAYS
ERxx = ER2(2:2:Nx2, 1:2:Ny2);
ERyy = ER2(1:2:Nx2, 2:2:Ny2);
URzz = UR2(2:2:Nx2, 2:2:Ny2);
SIGxx = SIG2(2:2:Nx2, 1:2:Ny2);
SIGyy = SIG2(1:2:Nx2, 2:2:Ny2);



% SHOW DEVICE
subplot(1,5,1);
imagesc(xa, ya, ER2.');
set(gca,'fontsize',20)
axis equal tight
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
ursrc = URzz(1, ny_ref);
urtrn = URzz(1, ny_trn);
ersrc = ERxx(1, ny_ref);
ertrn = ERxx(1, ny_trn);
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
ETAsrc = sqrt(URzz(1, ny_src)/ERxx(1, ny_src));   % Relative impedance of the source
nsrc = sqrt(URzz(1, ny_src)*ERxx(1, ny_src));     % Relative refrative index of the source
delt = 0.5*dt + 0.5*nsrc*dy/c0;                   % Times between the magnetic and electric field
t = [0:STEPS-1]*dt;
Hzsrc = exp(-((t-t0)/tau).^2);
Exsrc = -ETAsrc * exp(-((t-t0- delt)/tau).^2);

% COMPENSATE FOR DISPERSION 
k0 = 2*pi*f0/c0;                % Free sapce wave number
f = c0*dt/sin(pi*f0*dt)*sqrt((sin(k0*dy/2)/dy)^2);            
URzz = f*URzz;
ERxx = f*ERxx;
ERyy = f*ERyy;

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

% CALCULATE UPDATE COEFFICIENTS FOR Dx
sigy = sigy2(2: 2: Nx2, 1: 2: Ny2);
bDxy = exp(-(sigy/pml_ky + pml_ay)*dt/e0);
cDxy = sigy./(sigy*pml_ky + pml_ay*pml_ky^2).*(bDxy - 1);


% CALCULATE UPDATE COEFFICIENTS FOR Bz
sigy = sigy2(2: 2: Nx2, 2: 2: Ny2);
bBzy =  exp(-(sigy/pml_ky + pml_ay)*dt/e0);
cBzy = sigy./(sigy*pml_ky + pml_ay*pml_ky^2).*(bBzy - 1);

% CALCULATE UPDATE COEFFICIENTS FOR Ex and Ey
Ax = (dt*c0)*SIGxx;
Ay = (dt*c0)*SIGyy;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALISE FDTD TERMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALISE FILEDS TO ZERO 
Bz = zeros(Nx, Ny);
Hz = zeros(Nx, Ny);
Dx = zeros(Nx, Ny);
Dy = zeros(Nx, Ny);
Ex = zeros(Nx, Ny);
Ey = zeros(Nx, Ny);

% INITIALISE DERIVATIVE ARRAYS
dHzx = zeros(Nx, Ny);
dHzy = zeros(Nx, Ny);
dExy = zeros(Nx, Ny);
dEyx = zeros(Nx, Ny);

% INITIALISE CONVOLUTIONS
psiDx_ylo = zeros(Nx, NPML(1));
psiDx_yhi = zeros(Nx, NPML(2));
psiBz_ylo = zeros(Nx, NPML(1));
psiBz_yhi = zeros(Nx, NPML(2));

% INITIALISE INTEGRATION
IEx = zeros(Nx,Ny);
IEy = zeros(Nx,Ny);

% INITIALISE FOURIER TRANSFORM
K = exp(-1i*2*pi*FREQ*dt);                      % Calculate the kernals for each position on grid
HREF = zeros(Nx, NFREQ);                        % Reflected steady state fields FT
HTRN = zeros(Nx, NFREQ);                        % Transmitted steady state fields FT
SRC = zeros(1, NFREQ);                          % FT of the source

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN FDTD LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% MAIN LOOP -- ITERATE OVER TIME
%
for T = 1 : STEPS

    % CALCULATE DERIVATIVES OF E
        % dEyx
        for nx = 1 : Ny
            for nx = 1 : Nx - 1 
                dEyx(nx, ny) = (Ey(nx+1, ny) - Ey(nx, ny))/dx;
            end
            dEyx(Nx, ny) = (Ey(1, ny) - Ey(Nx, ny))/dx;
        end
        % dExy
        for nx = 1 : Nx
            for ny = 1 : Ny - 1 
                dExy(nx, ny) = (Ex(nx, ny+1) - Ex(nx, ny))/dy;
            end
            dExy(nx, Ny) = (Ex(nx, 1) - Ex(nx, Ny))/dy;
        end
    % INCOPORATE TF/SF CORRECTION
    dExy(:, ny_src-1) = dExy(:, ny_src-1) - Exsrc(T)/dy;
    
    % UPDATE CONVOLUTIONS FOR B FIELD UPDATES
    psiBz_ylo = bBzy(:, 1: NPML(1)).*psiBz_ylo ...
        + cBzy(:, 1: NPML(1)).*dExy(:, 1: NPML(1));
    psiBz_yhi = bBzy(:, Ny - NPML(2) + 1: Ny).*psiBz_yhi ...
        + cBzy(:, Ny - NPML(2) + 1: Ny).*dExy(:, Ny - NPML(2) + 1: Ny);
    
    % UPDATE B FROM E
    Bz = Bz - M*(dEyx - dExy/pml_ky);
    Bz(:, 1: NPML(1)) = Bz(:, 1: NPML(1)) ...
        + M*psiBz_ylo;
    Bz(:, Ny - NPML(2) + 1: Ny) = Bz(:, Ny - NPML(2) + 1: Ny) ...
        + M* psiBz_yhi;

    % UPDATE H from B
    Hz = Bz./URzz;

    % CALCULATE DERIVATIVES OF H
        % dHzx
        for ny = 1 : Ny
            dHzx(1, ny) = (Hz(1, ny) - Hz(Nx, ny))/dx;
            for nx = 2 : Nx
                dHzx(nx, ny) = (Hz(nx, ny) - Hz(nx-1, ny))/dx;
            end
        end
        % dHzy
        for nx = 1 : Nx
            dHzy(nx, 1) = (Hz(nx, 1) - Hz(nx, Ny))/dy;
            for ny = 2 : Ny 
                dHzy(nx, ny) = (Hz(nx, ny) - Hz(nx, ny-1))/dy;
            end
        end
    
    % INCOPORATE TF/SF CORRECTION
    dHzy(:, ny_src) = dHzy(:, ny_src) - Hzsrc(T)/dy;


    % UPDATE CONVOLUTIONS FOR D FIELD UPDATES
    psiDx_ylo = bDxy(:, 1: NPML(1)).*psiDx_ylo ...
       + cDxy(:, 1: NPML(1)).*dHzy(:, 1: NPML(1));
    psiDx_yhi = bDxy(:, Ny - NPML(2) + 1: Ny).*psiDx_yhi ...
       + cDxy(:, Ny - NPML(2) + 1: Ny).*dHzy(:, Ny - NPML(2) + 1: Ny);
    
    % UPDATE D FROM H
    Dx = Dx + M*(dHzy/pml_ky);
    Dx(:, 1: NPML(1)) = Dx(:, 1: NPML(1)) + M*psiDx_ylo;
    Dx(:, Ny - NPML(2) + 1: Ny) = Dx(:, Ny - NPML(2) + 1: Ny)...
        + M*psiDx_yhi;
    Dy = Dy - M*dHzx;
    % UPDATE E from D
    IEx = IEx + Ex;
    IEy = IEy + Ey;
    Ex = (Dx - Ax.*IEx)./ERxx;
    Ey = (Dy - Ay.*IEy)./ERyy;

    % UPDATE FOURIER TRANSFORM
    for nfreq = 1 : NFREQ
        HREF(:,nfreq) = HREF(:,nfreq) + (K(nfreq)^T)*Hz(:, ny_ref);
        HTRN(:,nfreq) = HTRN(:,nfreq) + (K(nfreq)^T)*Hz(:, ny_trn);
        SRC(nfreq) = SRC(nfreq) + (K(nfreq)^T)*Hzsrc(T);
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
            href = HREF(:,nfreq)/SRC(nfreq);
            aref = fftshift(fft(href))/Nx;
            RDE = abs(aref).^2.*real(kyref/kyinc);
            REF(nfreq) = sum(RDE); 

             % TRANSMISSION 
            htrn = HTRN(:,nfreq)/SRC(nfreq);
            atrn = fftshift(fft(htrn))/Nx;
            TDE = abs(atrn).^2.*real(ersrc/ertrn*kytrn/kyinc);
            TRN(nfreq) = sum(TDE);
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
        imagesc(xa, ya, Hz.');
        set(gca,'fontsize',20)
        axis equal tight;
        colorbar;
        caxis(d*[-1, 1]);  

        % SHOW FREQUENCY RESPONSE 
        subplot(1, 5, 3:5)
        plot(FREQ*10^(-12), CON, ':k', 'LineWidth',3);
        hold on;
        plot(FREQ*10^(-12), REF, '-r', 'LineWidth',3);
        plot(FREQ*10^(-12), TRN, '-b', 'LineWidth',3);
        plot(FREQ*10^(-12), ABSOR, '-G', 'LineWidth',3);
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
