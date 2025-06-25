clear; close all; clc;
%% ----- Globally defined parameters ---------
global h Lx Ly Nx Ny Nb ds s a c xu xv yu yv vS vN epsilon Uint Vint Xel Xer
global rhof rhos rhod mu Kb  alpha beta Kend NyPart CbFactor rhoMin
global umax L0 E hstar  Xib Xiu sigmau sigmab rho0 eta Cb rho Re gamma l0
global alpha0 beta0 Cb0 rho_0init Xib0 Xiu0 w0 theta lx0

%% ---------Function for angular displacement, velocity and acceleration ------------
fdisp = @(t,A,tau, k) ((4*A*(tau/k))./(1 + exp(-k*(t/tau - 0.5))));
fvel = @(t,A,tau,k)( 4*A*exp(-k*(t/tau - 0.5))./(1 + exp(-k*(t/tau - 0.5)) ).^2 );
facc = @(t,A,tau,k)((4*A*(k/tau))*(exp(-2*k*(t/tau-0.5))-exp(-k*(t/tau-0.5)))./(1+exp(-k*(t/tau-0.5))).^3);
%%%%%%%%%%%%======================================================================
%% ------------ Parameter values --------------------
NyPart = 18;                  %% Partition along the y axis

CbFactor = 100;               %% Magnification of adhesion bonds Cb


Ly0 = 2e-3;                   %% Height of fluid above eyewall (in m)
Lx0 = 3*Ly0;                  %% Half the domain size along x (in m) 
L0 = Ly0;                     %% Characteristics length
Lx = Lx0/L0;                  %% Dimensionless Length of half domain along x 
Ly = Ly0/L0;                  %% Dimensionless Length of domain along y  

k = 15;                       %% Parameter controling steepness of angular displacement
R = 1.2e-2;                   %% Radius of the eyeball
tau = [0.04, 0.04];           %% Periods of saccadic eye rotations
A =  [750,-750];              %% Peak angular velocities for various saccades
dt = 2e-4/NyPart;             %% Actual Time step 
tf = sum(tau);                %% Actual simulation time in seconds 
t = 0:dt:tf;                  %% Actual Time parttitions
Const_p = 1;                  %% Constant pressure intruoduce to account for Neumann BC
Num_SaveData = 100;           %% Number of data saved for the GIF plot
Nt = 50;                      %% Number of iteration for the p1_scheme

LocPar = [dt,tf, Const_p, Num_SaveData,Nt]; 

%%%%%%%%%%%%======================================================================
%% -------Computing eyewall displacement, velocity and accelation ------------------
Num_Saccade = length(A);        %% Number of Saccadic eye rotations
Theta = zeros(length(t),1);     %% Initiallized vector for angular displacement
Omega = zeros(length(t),1);     %% Initiallized vector for angular velocity
Alpha = zeros(length(t),1);     %% Initiallized vector for angular accleraction

Theta0 = 0;                     %% Initial displacement
start_time = 0;                 %% starting time
for n = 1:Num_Saccade
    % Determine the time range for the current segment
    end_time = start_time + tau(n);
    idx = ( t > start_time) & (t <= end_time);    %% Logical index for current segment

    %%% Compute angular displacement (theta), velocity (omega) and aacceration(alpha) for this segment
    %%% Note: We need to shift time to start at 0 for each segment
    Theta(idx) = fdisp(t(idx) - start_time, A(n), tau(n), k) + Theta0;  %% Angular displacement

    Omega(idx) = fvel(t(idx) - start_time, A(n), tau(n), k);            %% Angular Velocity

    Alpha(idx) = facc(t(idx) - start_time, A(n), tau(n), k);            %% Angular Acceleration

    %%% Update the initial displacement for the next segment 
    Theta0 = Theta(find(idx, 1, 'last'));  

    %%% Update the start time for the next segment
    start_time = end_time;
end

%%% ----linear displacement, velocity and accelration------------ 
S = Theta*(pi/180)*R;           %% Linear displacement
WallVel = Omega*(pi/180)*R;     %% Linear velocity
WallAcc = Alpha*(pi/180)*R;     %% Linear acceleration

%%% ----Nondimensionalizing linear displacement, velocity and accelration------------
Vmax = max(abs(WallVel));        %% maximum linear velocity
amax = max(abs(WallAcc));        %% maximum linear acceleration 
S = S/L0;                        %% Dimensionless Lineae displacement 
WallVel = WallVel/Vmax;          %% Dimensionless Linear velocity
WallAcc = WallAcc/amax;          %% Dimensionless Linear acceleration

%%%%%%%%%%%%======================================================================
%%%%%------- These model parameters are kept fixed -----------------------------------
hstar = 2.0e-4;                %% tickness of retina in [m]
umax = Vmax;                   %% maximum wall velocity in  [m/sec]       
E = 1.21e3;                    %% Yong mudulus in [N/m^2]; 
mu = 1.065e-3;                 %% Dynamic viscosity in [Ns/m^2]
rhof = 1000;                   %% Fluid density in [kg/m^3]
rhos = 1040;                   %% material/retina density in [kg/m^3]
rhod = (rhos - rhof)*hstar;    %% Difference between material and fluid density
Kb = (E*hstar^3)/12;           %% Bending stiffness in [Nm^2]

%%%%%--These model parameters can be varied to until desire result is achieved-----------
alpha0 = -100;                 %% Coefficient of Goldstein feedback law  in [N/m^4]      
beta0 = -10;                   %% Coefficient of Goldstein feedback law  in [Ns/m^4]
Cb0 = 2.0;                     %% strength of adhesion bond in N/m                 
rho_0init = 100;               %% Density of available linker proteins [1/m]
Kb = (E*hstar^3)/12;           %% Bending stiffness

Xib0 = 600;                    %% Binding rate   
Xiu0 =   0.6;                  %% Unbinding rate
sigmab = 5e5*L0;               %% exponent of binding term
sigmau = 1.8e5*L0;             %% 2e5*L0 exponent of binding term  

%%%%%--These model parameters can be varied to until desire result is achieved----------- 
alpha = (alpha0*L0^2)/(rhod*umax^2);     %% Dimensioless G. Feedback coefficient    
beta = (beta0*L0)/(rhod*umax);           %% Dimensioless G. Feedback coefficient
Cb = (Cb0*L0)/(rhod*umax^2);             %% Dimensioless adhession bond strength
gamma = (Kb)/(rhod*umax^2*L0^2);         %% dimensionless bending stiffness
Re = (rhof*umax*L0)/(mu);                %% Reynolds number
rho = rhod/(rhof*L0);                    %% dimensionless coefficient of fluid force
rho0 = rho_0init*L0;                     %% Dimensioless available linker protein
Xiu = Xiu0*L0/umax;                      %% Dimensionless unbinding rate
Xib =  Xib0*L0/umax;                     %% Dimensionless binding rate           
Kend = gamma;                            %% spring constant connecting the detached ends

%%% -------------- Space parameters (dimensionless) ---------------------
h =  Ly/NyPart;               %% spatial resolution                 
Nx = (2*Lx)/h;                %% Number of grid cells along x 
Ny = Ly/h;                    %% Number of grid cells along y
Uint = Ny*(Nx+1);             %% Size of interior nodes for u
Vint = (Ny-1)*(Nx);           %% Size of interior nodes for v

%%% -------------- Initial condition and other parameters ---------------------
theta = (pi/180)*7;                   %% angle of elevation of the bulge
l0 = 0.70e-3/L0;                      %% initial detached length
w0 = Ly/100;                          %% natural length of the linker proteins
c = l0*sin(theta);                    %% heigh of the bulge above the eyewall
epsilon =2e-2 + l0*(1 - cos(theta)) ; %% radius of retinal hole
lx0 = l0*cos(theta);                  %% x-component of the initial detached length
a = lx0/(sqrt(-log(w0/c)));           %% Exponential constant for the initial detached retina
rhoMin = rho0/100;                    %% Maximum displacement beyond which detachment occur

eta = 2e4*L0;                         %% steepness of rhob at the clamp end [dimensionless]

vS = 0;                                %% Boundary condition for v at the bottom                     
vN = 0;                                %% Boundary condition for v at the top  
ds = h;                                %% step size for Lagrangian coordinates 
s = epsilon : ds : Lx;                 %% Lagrangian partition for half side
Nb = length(s);                        %% Number of Lagrangian nodes 
xu = -Lx : h : Lx;                     %% x coordinate for the u verlocity
xv = -Lx-h/2 : h : Lx+h/2;             %% x coordinate for the v verlocity
yu = -h/2: h: Ly+h/2;                  %% y coordinate for the u verlocity
yv = 0 : h : Ly;                       %% y coordinate for the v verlocity

%% -------- Defining initial velocity and Lagrangian coordiate -----------
u0 = zeros(Ny+2, Nx+1);                %% Initial u-velocity
u0(1,:) = 2*WallVel(1) - u0(2,:);      %% specifying lower boundary velocity for u 
v0 = zeros(Ny+1, Nx+2);                %% Initial u-velocity
[Xel, Xer] = LagStructure();           %% Initial geometry of the detached retina

%%%%%%%%%%============= OUTPUT===================================
[Xl, Xr] = LagStructure();             %% Initial geometry of the detached retina

tic
[pc0,pc1, pp,pp1,pd1,pd2, pc, pdet, pu, pv]=MainFlow(LocPar,u0,v0,S,Theta,  WallAcc, WallVel,Omega);
Time = toc

% S = S*L0;
% WallVel = WallVel*umax; 
% 
% plot(t, Theta, 'b', 'LineWidth',1.5)