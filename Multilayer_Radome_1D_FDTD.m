% Multilayer_Radome_1D_FDTD.m
% Calculate Reflectance and Transmittance of a multilayer radome

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
meters = 1;
seconds = 1;

centimeters  = 1e-2 * meters;
millimeters = 1e-3 * meters;
inches      = 2.54 * centimeters;
hertz       = 1/seconds;
kilohertz   = 1e3 * hertz;
megahertz   = 1e6 * hertz;
gigahertz   = 1e9 * hertz;

% CONSTANTS
c0 = 299792458 * meters/seconds;
e0 = 8.8541878176e-12 * 1/meters;
u0 = 1.2566370614e-6 * 1/meters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1 -- DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE PARAMETERS
fmin = 17.0 * gigahertz;
fmax = 32.0 * gigahertz;
NFREQ = 500;
FREQ = linspace(fmin,fmax,NFREQ);

% DEVICE PARAMETERS FOR MULTILAYER RADOME
er1 = 1.1; % Millifoam
er2 = 4.0; % GFRP
er3 = 1.0; % Air
L1  = 2.1 * millimeters;
L2  = 0.24 * millimeters;
L3  = 0.48 * millimeters;

urmax = 1;

% GRID PARAMETERS
ermax = max([er1 er2 er3]);    % Find greatest dielectric constant
nmax  = sqrt(ermax * urmax);   % Greatest dk used to calculate refractive index

DDAT  = [L2 L1 L3 L1 L2];      % Setup material length array 
erDAT = [er2 er1 er2 er1 er2]; % Setup material dielectric constant array

NRES_LAM = 60;                 % Resolution in cells per wavelength
NRES_D   =2;                   % ??
NSPC     = [100 100];          % Buffer of free space before and after device

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2 -- CALCULATE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE INITIAL GRID RESOLUTION
Lmin    = min([L1 L2 L3]);  % Find minimum physical length
lam_min = c0/fmax/nmax;     % Calculate minimum lambda at fmax and in nmax
dz1     = lam_min/NRES_LAM; % Minimum dz based on resolution min lam / resolution
dz2     = Lmin/NRES_D;      % Find min dz based on Lmin physical dimention
dz      = min([dz1 dz2]);   % Choose dz to be smallest of dz1 and dz2

% SNAP GRID TO CRITICAL DIMENSIONS
nz = ceil(Lmin/dz)
dz = Lmin/nz

% CALCULATE NUMBER OF GRID CELLS
Nz = ceil(L2/dz) + ceil(L1/dz) + ceil(L3/dz) +ceil(L1/dz) + ceil(L2/dz)+ sum(NSPC) + 3

% CALCULATE ARRAY OF E FIELD POSITION
za = [0:Nz-1]*dz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3 -- MODEL DEVICE ON THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATERIALS
ERxx = er3*ones(1,Nz);  % er3 is free-space permitivity
URyy = ones(1,Nz);      % ur =  1 since non magnetic material

% CALCULATE START AND STOP INDICIES OF SLAB & INCORPORATE SLAB INTO ERxx
nz0 = 2 + NSPC(1) + 1;           % nz0 is start of radome
nz1 = nz0;                       % set nz1 to start index of radome
for nd = 1 : length(DDAT)
    nz2 = nz0 + round(sum(DDAT(1:nd))/dz) - 1; % Loop for start/stop of each layer in DDAT
    ERxx(nz1:nz2) = erDAT(nd);        % Insert er for material
    nz1 = nz2 + 1;
end

 plot(za/millimeters,ERxx) % Uncomment to plot er of structure layers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4 -- CALCULATE THE TIME STEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DETERMINE REFRACTIVE INDEX AT BOUNDARIES
nbc = sqrt(URyy(1)*ERxx(1));

% CALCULATE TIME STEP FOR PABC
dt = nbc*dz/(2*c0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5 -- CALCULATE THE SOURCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% POSITION OF SOURCE
k_src = 2;    % 2nd cell into simulation domain, from the left.

% CALCULATE GAUSSIAN PULSE PARAMETERS
tau = 0.5/fmax;
t0 = 3*tau;

% CALCULATE TIME ARRAY
N = 2000;
t = [0:N-1]*dt;

% CALCULATE Ex SOURCE
Exsrc = exp(-((t - t0)/tau).^2);

% CALCULATE Hy SOURCE
nsrc  = sqrt(ERxx(k_src)*URyy(k_src));
s     = 0.5*nsrc*dz/c0 - dt/2;
A     = sqrt(ERxx(k_src)/URyy(k_src));
Hysrc = A*exp(-((t - t0 + s)/tau).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 6 -- CALCULATE UPDATE COEFFICIENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE UPDATE COEFFICIENTS
mEx = -(c0*dt/dz)./ERxx;
mHy = -(c0*dt/dz)./URyy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 7 -- INITIALIZE FOURIER TRANSFORMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE KERNALS
K = exp((-1i*2*pi*dt)*FREQ);

% INITIALIZE FOURIER TRANSFORM ARRAYS
ExR = zeros(1,NFREQ);
ExT = zeros(1,NFREQ);
SRC = zeros(1,NFREQ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 8 -- INITIALIZE FIELDS TO ZERO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE FIELDS TO ZERO
Ex = zeros(1,Nz);
Hy = zeros(1,Nz);

% INITIALIZE BOUNDARY FIELD TERMS
E1 = 0;    E2 = 0;
H1 = 0;    H2 = 0;

%
% MAIN FDTD LOOP -- ITERATE OVER TIME
%
for n = 1 : N
  
  % Step a -- Update Ex Boundary Terms
  E2 = E1;
  E1 = Ex(Nz);
  
  % Step b -- Update Ex
  Ex(1) = Ex(1) + mEx(1)*(Hy(1) - H2);  % With Direchlet BC so Hy(0) = 0
  for k = 2 : Nz
    Ex(k) = Ex(k) + mEx(k)*(Hy(k) - Hy(k-1));
  end
  
  % Step c -- Incorporate TF/SF Correction Term
  Ex(k_src) = Ex(k_src) - mEx(k_src)*Hysrc(n);
  
  % Step d -- Update Hy Boundary Terms
  H2 = H1;
  H1 = Hy(1);
  
  % Step e -- Update Hy
  
  for k = 1 : Nz -1
    Hy(k) = Hy(k) + mHy(k)*(Ex(k+1) - Ex(k));
  end
  Hy(Nz) = Hy(Nz) + mHy(Nz)*(E2 - Ex(Nz));  % With Direchlet BC so Ex(k+1) = 0
  
  % Step f -- Incorporate TF/SF Correction Term
  Hy(k_src-1) = Hy(k_src - 1) - mHy(k_src - 1)*Exsrc(n);
  
  % Step g -- Update Fourier Transforms
  for nf = 1 : NFREQ
    ExR(nf) = ExR(nf) + (K(nf)^n)*Ex(1);
    ExT(nf) = ExT(nf) + (K(nf)^n)*Ex(Nz);
    SRC(nf) = SRC(nf) + (K(nf)^n)*Exsrc(n);  
  end
  
  % Step h -- Visualize Simulation
  if ~mod(n,50)
    
    % Calculate Spectra
    REF = abs(ExR./SRC).^2;
    TRN = abs(ExT./SRC).^2;
    CON = REF + TRN;
    
    % Prepare Figure Window
    clf;
    % Show Fields
    subplot(212);
    plot(za/millimeters,Ex, '-b');
    hold on;
    plot(za/millimeters,Hy, '-r');
    hold off;
    xlim([za(1) za(Nz)]/millimeters/meters);
    ylim([-1.8 1.8]);
    xlabel('z-axis');
    title(['Step ' num2str(n) ' of ' num2str(N) ]);
    % Show Spectra
    subplot(211);
##    plot(FREQ/gigahertz,REF, '-r');
##    hold on;
##    plot(FREQ/gigahertz,TRN, '-b');
##    plot(FREQ/gigahertz,CON, ':k');
plot(FREQ/gigahertz,10*log10(REF), '-r');
hold on;
plot(FREQ/gigahertz,10*log10(TRN), '-b');
plot(FREQ/gigahertz,10*log10(CON), ':k');
    
    hold off;
    xlim([fmin fmax]/gigahertz);
    ylim([-80 1.3]); % Change limits to suit...
    xlabel('Frequency (GHz)');
    title('Spectra');
    drawnow;
  end
  
end

