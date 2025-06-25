function [pc0,pc1, pp,pp1,pd1,pd2, pc, pdet, pu, pv]=MainFlow(LocPar,u0,v0,S,Theta, WallAcc,WallVel,Omega)

global h Lx Ly Nx Ny Nb ds s a c xu xv yu yv vS vN epsilon Uint Vint Xel Xer
global rhof rhos rhod mu Kb  alpha beta Kend NyPart CbFactor rhoMin
global umax L0 E hstar  Xib Xiu sigma rho0 eta Cb rho Re gamma l0
global alpha0 beta0 Cb0 rho_0init Xib0 Xiu0 w0 ChibFactor

%%% ------ Local parameters ---------------
% LocPar = [dt,tf, Const_p, Nt, Num_SaveData];
dt = LocPar(1);               %% dimensional time step
tf = LocPar(2);               %% dimensional final simulation time
Const_p = LocPar(3);          %% Constant pressure intruoduce to account for Neumann BC
Num_SaveData = LocPar(4);     %% Number of saved data for GIF plot
Nt = LocPar(5);               %% Number of iteration for the explicity 

dt = dt*umax/L0;              %% Dimensionless time step
tf = tf*umax/L0;              %% Dimensionless simulation time
time = 0:dt:tf;               %% Dimensionless time

NumOfIter = floor(tf/dt);                   %% maximum number of time iterations
freqSave = floor(NumOfIter/Num_SaveData);   %% How often the data is saved

%%% ****** Simulation to get the second intial condition ******************
dt1 = dt/(Nt);                      %% time step for the explicit scheme
[un, vn, Xl, Xr,p0,rhoL0, rhoR0]=P1FirstIter(u0,v0,dt1,Const_p, WallAcc,WallVel,Nt);

%%%%%%---------User define functions ----------------------------------
[Xl0, Xr0] = LagStructure();
[Ap] = PLapMat(Const_p);
[Au,Av] = UV_LapMat(dt);
[FaL0, FaR0] = Falpha0(u0, v0, Xl, Xl0, Xr, Xr0,dt);
[wL0, wR0] = BondDispW(Xl, Xr);

%%% Computing the initial detached length 
[wLinit, wRinit] = BondDispW(Xl0, Xr0);
[RhobLinit, RhobRinit] = InitRhob(Xl0, Xr0);     %% rho_b at time zero
[distL0, distR0] = DetachRate(wLinit, wRinit,RhobLinit, RhobRinit, Xl0, Xel, Xr0, Xer, dt);


%%%% ----------------------------------------------------------------------------
SaveMag = zeros(Ny+1, Nx+1, Num_SaveData + 1);    %% the magnite of velocity
SaveU = zeros(Ny+1, Nx+1, Num_SaveData + 1);      %% u-component of the velocity
SaveV = zeros(Ny+1, Nx+1, Num_SaveData + 1);      %% v-component of the velocity
Vort = zeros(Ny+1, Nx+1, Num_SaveData + 1);       %% the vorticity of the flow
SaveXL = zeros(Nb, 2, Num_SaveData+1);            %% Left Lagrangian structure
SaveXR = zeros(Nb, 2, Num_SaveData+1);            %% Right the Lagrangian structure
SaveRho = zeros(Nb, 2, Num_SaveData+1);           %% Density of linker proteins 
Converge = zeros(NumOfIter-2, 2);                 %% Continuity data
DLent = zeros(NumOfIter-1, 3);                    %% Data for detached length
CFL_Cond = zeros(NumOfIter-2, 3);                 %% CFL condition
ClampT = zeros(NumOfIter-2, 3);                   %% Total force at the clamp edge
ClampP = zeros(NumOfIter-2, 3);                   %% Lagrangian force at the clamp edge
SnapTime = zeros(Num_SaveData+1,1);                    %% Selected time points for clock

Udata0 = 0.5*( u0(2:end,:) + u0(1:end-1, :) );
Vdata0 = 0.5*( v0(:,2:end) + v0(:, 1:end-1) );
magU0 = sqrt(Udata0.^2+Vdata0.^2 + eps);
vort0 = (v0(:, 2:end) - v0(:,1:end-1) )/h - ( u0(2:end,:) - u0(1:end-1,:) )/h;

SaveMag(:,:,1) = magU0;
SaveU(:,:,1) = Udata0;
SaveV(:,:,1) = Vdata0;
Vort(:,:,1) = vort0;
SaveXL(:,:,1) = Xel;
SaveXR(:,:,1) = Xer;
SaveRho(:, :, 1) = [RhobLinit, RhobRinit];
DLent(1,:) = [time(1), distL0, distR0];


count = 2;

for step = 3: NumOfIter

    ustar = zeros(Ny+2, Nx+1);
    vstar = zeros(Ny+1, Nx+2);
    pprime = zeros(Ny+2, Nx+2);
    Grad2P = zeros(Ny+2, Nx+2);   %% for computing pressure at half time step

    uBC = WallVel(step);          %% Dirichlet boundary condition for the fluid velocity

    %% STEP 1: -- Solving for intermediate velocity u* and v*------------
    [RHSu, RHSv]= P2_RHS(FaL0, FaR0,u0,un,v0,vn,p0,Xl,Xl0, Xr,Xr0,dt, uBC);

    %%% ---Solving the Poisson equation for u* and applying the BC -----
    us_int = Au\RHSu;                                       %% Computing interior values for u*
    ustar(2:end-1, :) = reshape(us_int, Nx+1, Ny)';         %% Converting interior values for u* into matrix
    ustar(1,:) = 2*uBC - ustar(2,:);                         %% Dirichlet BC at the buttom
    ustar(end,:) = ustar(end-1,:);                          %% Neumann BC at the top
    ustar(:,1) = ustar(:,end);   ustar(:,end) = ustar(:,1); %% Periodic BC on the left and right

    %%% ---Solving the Poisson equation for v* and applying the BC ------
    vs_int = Av\RHSv;                                        %% Computing interior values for v*
    vstar(2:end-1, 2:end-1) = reshape(vs_int, Nx,Ny-1)';     %% Converting interior values for u* into matrix
    vstar(:,1) = vstar(:,end-1);              %% periodic BC on the left; vstar(:,1) = vstar(:,end-1)
    vstar(:,end) = vstar(:,2);                %% periodic BC on the right; vstar(:,end) = vstar(:,2)
    vstar(1,:) = vS;    vstar(end,:) = vN;    %% Dirichlet BC at the top and Bottom

    %% --- STEP 2: Computing the pressure correction p' -----
    %%% -- Here, we solve grad^2(p') = (rho/dt)*Div(U*) where U* = (u*, v*)--
    dustar_dx = ( ustar(:, 2:end) - ustar(:, 1:end-1) )/h;  %% derivative of u* wrt x
    dvstar_dy = ( vstar(2:end, :) - vstar(1:end-1 , :) )/h; %% derivative of v* wrt y
    gradU = (dustar_dx(2:end-1, :) + dvstar_dy(:, 2:end-1)); %% u*_x+v*_y for interior
    gradU_vec = reshape(gradU', (Ny)*(Nx),1);                %% shaping into a vector

    pint = Ap\((1/dt)*(gradU_vec));            %% pressure correction at the interior nodes
    pprime(2:end-1, 2:end-1) = reshape(pint, Nx, Ny)';
    pprime(1,:) = pprime(2,:);                 %% applying left BC
    pprime(end,:) = pprime(end-1,:);           %% applying right BC
    pprime(:,1) = pprime(:,2);                 %% applying down BC
    pprime(:, end) = pprime(:,end-1);           %% applying top BC

    %% -- STEP 3: Finding the velocity at the next time level ----
    dpprime_dx = (pprime(:,2:end) - pprime(:,1:end-1))/h;
    dpprime_dy = ( pprime(2:end,:) - pprime(1:end-1, :) )/h;
    %%% Computing u & v velocities at time level n+1
    unew = ustar -dt*dpprime_dx;     %% u velocity at time n+1
    unew(1,:) =  2*uBC - unew(2,:);                        %% Applying down BC
    unew(end,:) = unew(end-1,:);
    unew(:,1) = unew(:,end);   unew(:,end) = unew(:,1);

    vnew = vstar - dt*dpprime_dy;         %% v velocity at time n+1
    vnew(:,1) = vnew(:,end-1);                  %% enforcing the periodic BC on the left
    vnew(:,end) = vnew(:,2);                    %% enforcing the periodic BC on the right
    vnew(1,:) = vS;    vnew(end,:) = vN;

    %%% ---Computing gradient(U) to check whether the continuity equation is satisfied or not------
    Ux = (unew(:, 2:end) - unew(:, 1:end-1))/(h);
    Vy = (vnew(2:end,:) - vnew(1:end-1, :) )/(h);
    divU = Ux(2:end-1, :) + Vy(:, 2:end-1);
    Converge(step-2,:) = [time(step), norm(divU,"fro")];

    %% -- STEP 4: Computing the pressure p^(n+0.5) ----
    %%% NOTE: P^(n+0.5) = P^(n-0.5) + Pprime - (dt)/(2*Re)*Grad^2(Pprime)
    p_xx = (pprime(:, 3:end) - 2*pprime(:, 2:end-1) + pprime(:,1:end-2))/(h^2);
    p_yy = (pprime(3:end, :) - 2*pprime(2:end-1, :) + pprime(1:end-2,:))/(h^2);

    Grad2P(:,2:end-1) = p_xx;
    Grad2P(2:end-1,:) = Grad2P(2:end-1,:) + p_yy;
    p_half = p0 + pprime - (dt/(2*Re))*Grad2P;       %% pressure at next half time step

    p_half(:,1) = p_half(:,2);
    p_half(:,end) = p_half(:,end-1);
    p_half(1,:) = p_half(2,:);
    p_half(end,:) = p_half(end-1, :);

    %% --STEP 5: Extrapolating the pressure at next time level n+1 ---
    %%% NOTE: P^(n+1) = 3/2*P^(n+0.5) - 1/2*P^(n-0.5)
    pnew = (3/2)*p_half - (1/2)*p0;

    %%%% ----- Computing new Lagrangian position --------
    [XLnew, XRnew,cfl] = UpdateLag(WallAcc(step), wL0, wR0, rhoL0, rhoR0, FaL0, FaR0, un,vn,Xl, Xl0, Xr,Xr0, dt);

    CFL_Cond(step-1, :) = [time(step), cfl];

    %% ------ Update initial data ----------------
    [UibL,UibR] = Uib(unew,vnew, XLnew, XRnew);

    [rhobL, rhobR] = RhoB(wL0, wR0, rhoL0, rhoR0, Xl, Xl0, Xr, Xr0, dt);
    %%% Computing the force at the clamp edge
    [CFtot, Lff] = ClampForce(WallAcc(step), wL0, wR0, rhobL,rhobR,rhoL0, rhoR0, FaL0, FaR0, un, vn, Xl,Xl0, Xr,Xr0,dt);
    ClampT(step-2, :) = [time(step), CFtot];
    ClampP(step-2, :) = [time(step), Lff];
    [distL, distR] = DetachRate(wL0, wR0,rhoL0, rhoR0, Xl, Xl0, Xr, Xr0, dt);
    DLent(step-1, :) = [time(step), distL, distR];

    dxLdt = (XLnew - Xl)/dt;
    dxRdt = (XRnew - Xr)/dt;

    FaL0 = FaL0 + alpha*(UibL - dxLdt )*dt;
    FaR0 = FaR0 + alpha*(UibR - dxRdt )*dt;
    [wL, wR] = BondDispW(XLnew, XRnew);

    rhoL0 = rhobL;
    rhoR0 = rhobR;

    u0 = un;        un = unew;
    v0 = vn;        vn = vnew;

    Xl0 = Xl;     Xl = XLnew;

    Xr0 = Xr;     Xr = XRnew;

    wL0 = wL;     wR0 = wR;

    p0 = p_half;

    if rem(step,freqSave) == 0
        udata = 0.5*( unew(2:end,:) + unew(1:end-1, :));
        vdata = 0.5*( vnew(:,2:end) + vnew(:,1:end-1));

        magU = sqrt(udata.^2 + vdata.^2 + eps);        %% compute the magnitude of the velocity

        SaveMag(:,:,count) = magU;
        SaveU(:,:,count) = udata;     %% %./magU;
        SaveV(:,:,count) = vdata;     %%./magU;
        Vort(:,:,count) =(vnew(:, 2:end)-vnew(:,1:end-1))/h - (unew(2:end,:) - unew(1:end-1,:))/h;
        SaveXL(:,:, count) = XLnew;
        SaveXR(:,:, count) = XRnew;
        SaveRho(:,:, count) =  [rhobL, rhobR];
        SnapTime(count) = time(step);

        count = count + 1;
    end
end

%%% Scaling the quantities back to the dimensional form 
SaveXL = SaveXL*L0;
SaveXR = SaveXR*L0;
SaveRho = SaveRho/L0;
SnapTime = SnapTime*(L0/umax);
SaveU = SaveU*umax;
SaveV = SaveV*umax;
Converge(:,1) = Converge(:,1)*(L0/umax);
DLent(:,1) =  DLent(:,1)*(L0/umax);
DLent(:,2:3) = DLent(:,2:3)*L0;
CFL_Cond(:,1) = CFL_Cond(:,1)*(L0/umax);
XLnew = XLnew*L0;
XRnew = XRnew*L0;
rhobL = rhobL/L0;
rhobR = rhobR/L0;
ClampP(:,1) = ClampP(:,1)*(L0/umax);
ClampT(:,1) = ClampT(:,1)*(L0/umax);
Xel = Xel*L0;
Xer = Xer*L0;
Lx = Lx*L0;
Ly = Ly*L0;
h = h*L0;
S = S*L0;
WallVel = WallVel*umax; 

%%%%%%%%%%==========================================================================
%% -------------- Visualizing the output -------------------------------------------
x = -Lx: h : Lx;                             %% x-values at the nodes
y = 0 : h : Ly;                              %% y-values at the nodes
[X,Y] = meshgrid(x,y);                       %% meshgrid for contour plot

%%%%%%%%**************************************************************************
%%%  --------GIF for the fluid flow-------------------
figure
filename1 = 'P2Flow.gif';
for k = 1: Num_SaveData+1
    pp = contourf(X,Y,SaveMag(:,:,k),"LineStyle","none");
    hold on
    plot(SaveXL(:,1,k), SaveXL(:,2,k), "r-",SaveXR(:,1,k), SaveXR(:,2,k), "r-","LineWidth",1.20);
    plot(Xel(:,1), Xel(:,2), "k-", Xer(:,1), Xer(:,2), "k-","LineWidth",1.20);
    colorbar

    sc = streamslice(X,Y, SaveU(:,:,k), SaveV(:,:,k), 6);
    set(sc,'LineWidth', 0.8,'color', 'w') %% [0.8, 0.55, 0.40]

    cmap = parula;                         % Get the parula colormap
    lighter_cmap = cmap*0.9 + 0.3;         % Blend with white (0.7 original + 0.3 white)
    lighter_cmap (lighter_cmap > 1) = 1;   % Keep values within valid range
    colormap(lighter_cmap);

    colorbar;

    grid on;
    xlim([-Lx,Lx]); ylim([0, Ly/2]);
    xlabel("\bf eyewall [m]"); ylabel("\bf vitreous [m]");
    title('\bf Flow profile');
    %%% Adding display time
    SimTime1 = sprintf('Time = %.2f s', SnapTime(k) );
    text(0.08*Lx, 0.25*Ly, SimTime1, 'FontSize',10, 'color','w','FontWeight', ...
        'bold','BackgroundColor','k' );
   
    pause(0.1);

    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,filename1,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename1,'gif','WriteMode','append');
    end
end

close(gcf);

%%%%%%%%**************************************************************************
%%%  ---------- GIF for the peeling of the retina -------------------
figure
PeelName2 = 'Detachment_Dynamics.gif';
for kk = 1: Num_SaveData+1
    pd1 = plot(Xel(:,1), Xel(:,2),'r', Xer(:,1), Xer(:,2), 'r','LineWidth',1.2);
    hold on
    pdel1 = plot(SaveXL(:,1,kk), SaveXL(:,2,kk),'b', SaveXR(:,1,kk), SaveXR(:,2,kk), 'b','LineWidth',1.2);

    title('\bf Detachment profile');
    grid on;
    xlim([-Lx,Lx]); ylim([0, Ly/5;]);
    xlabel("\bf eyewall [m]"); ylabel("\bf vitreous [m]");

    %%% Adding display time
    SimTime2 = sprintf('Time = %.2f s', SnapTime(kk) );
    text(0.08*Lx, 0.10*Ly, SimTime2, 'FontSize',10, 'color','w','FontWeight', ...
        'bold','BackgroundColor','k' );

    pause(0.1);

    drawnow
    frame1 = getframe(1);
    im = frame2im(frame1);
    [imind,cm] = rgb2ind(im,256);
    if kk == 1
        imwrite(imind,cm,PeelName2,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,PeelName2,'gif','WriteMode','append');
    end
    delete(pdel1);
end

close(gcf);

%%%%%%%%**************************************************************************
%%%  ----------GIF for the linker proteins-------------------
figure
DensName2 = 'Adhesion.gif';
for kr = 1: Num_SaveData+1
    pd2 = plot(Xel(:,1), SaveRho(:,1,1),'r', Xer(:,1), SaveRho(:,2,1), 'r','LineWidth',1.2);
    hold on
    pdel2=plot(SaveXL(:,1,kr), SaveRho(:,1,kr),'b', SaveXR(:,1,kr), SaveRho(:,2,kr), 'b','LineWidth',1.2);

    title('\bf Density of linker proteins');

    grid on;
    xlim([-Lx,Lx]);
    xlabel("\bf eyewall [m]"); ylabel("\bf \rho_b [m^{-1}]");
    %%% Adding display time
    SimTime3 = sprintf('Time = %.2f s', SnapTime(kr));
    text(0.3*Lx, 0.55*rho_0init, SimTime3, 'FontSize',10, 'color','w','FontWeight', ...
        'bold','BackgroundColor','k' );

    pause(0.1);

    drawnow
    frame1 = getframe(1);
    im = frame2im(frame1);
    [imind,cm] = rgb2ind(im,256);
    if kr == 1
        imwrite(imind,cm,DensName2,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,DensName2,'gif','WriteMode','append');
    end

    delete(pdel2);
end

close(gcf);

%%%%%%%%=======================================================================
%%% ------- Checking continuity of the velocity field-----------
figure
plot(Converge(:,1), Converge(:,2),"r", "LineWidth",1.0)
xlabel("\bf Time [s]");     ylabel("\bf continuity condtion, \nabla\cdot u")
pc= gcf;
exportgraphics(pc,'etdRdConvergence.png','Resolution',300);

%%%------ tracking the detached length over time ------------------
% Eftv_Left = DLent(end,2) - DLent(1,2);
% Eftv_Right = DLent(end,3) - DLent(1,3);

figure
plot(DLent(:,1), DLent(:,2), 'r--', DLent(:,1), DLent(:,3), 'b', 'LineWidth',2.0);
% title(['\bf Effective detach length [Left, Right] = ', num2str([Eftv_Left ,Eftv_Right])]);
xlabel('\bf Time [s]'); ylabel('\bf detached length [m]');
legend('left retina', 'right retina', 'location', 'southeast');
pdet= gcf;
exportgraphics(pdet,'detachedLength.png','Resolution',300);

%%% -----------Checking the CFL condition ------------------------
figure
semilogy(CFL_Cond(:,1), CFL_Cond(:,2), 'r', CFL_Cond(:,1), CFL_Cond(:,3), 'b', 'LineWidth',1.2);
% title('CFL condition for \rho_b over time');
xlabel('\bf Time [s]'); ylabel('\bf CFL value');
legend('left retina', 'right retina');
pc0= gcf;
exportgraphics(pc0,'CFL_Cond.png','Resolution',300);

%%% -----------Computing the total forces at the clamp edge ------------------------
figure
semilogy(ClampT(:,1), ClampT(:,2), 'r', ClampT(:,1), ClampT(:,3), 'b', 'LineWidth',1.2);
title('Total force at the clamp edge');
xlabel('\bf Time [s]'); ylabel('\bf clamp force [N]');
legend('left retina', 'right retina');
pc1= gcf;
exportgraphics(pc1,'Total_Clamp.png','Resolution',300);

%%% -----------Computing the Lagrangian forces at the clamp edge ------------------------
figure
semilogy(ClampP(:,1), ClampP(:,2), 'r', ClampP(:,1), ClampP(:,3), 'b', 'LineWidth',1.2);
title('Adhesion and fluid force at the clamp edge');
xlabel('\bf Time [s]'); ylabel('\bf clamp force [N]');
legend('left retina', 'right retina');
pcmp= gcf;
exportgraphics(pcmp,'Partial_Clamp.png','Resolution',300);

%%% ------- Final u and v velocity profile -----------
ufinal = 0.5*( unew(2:end,:) + unew(1:end-1,:));   %% final u velocity
vfinal = 0.5*( vnew(:,2:end) + vnew(:,1:end-1));   %% final u velocity

ufinal = ufinal*umax;
vfinal = vfinal*umax;
time = time*(L0/umax);

figure
contourf(X,Y,ufinal,51,"LineStyle","none");
hold on
plot(XLnew(:,1), XLnew(:,2), "m -",XRnew(:,1), XRnew(:,2), "m -","LineWidth",1.2);
plot(Xel(:,1), Xel(:,2), "g -", Xer(:,1), Xer(:,2), "g -","LineWidth",1.2);
title("final u-velocity profile");
xlabel("position");    ylabel("u-component ")
colormap('jet')
xlim([-Lx,Lx]), ylim([0,Ly/5]);
colorbar
pu= gcf;
exportgraphics(pu,'Final u.png','Resolution',300);
%%%%---------------v -velocity-----------------------------
figure
contourf(X,Y,vfinal,51,"LineStyle","none");
hold on
plot(XLnew(:,1), XLnew(:,2), "m -",XRnew(:,1), XRnew(:,2), "m -","LineWidth",1.2);
plot(Xel(:,1), Xel(:,2), "g -", Xer(:,1), Xer(:,2), "g -","LineWidth",1.2);
title("final v-velocity profile")
xlabel("position ")
ylabel("v-component")
colormap('jet')
xlim([-Lx,Lx]), ylim([0,Ly/5]);
colorbar
pv= gcf;
exportgraphics(pv,'Final v.png','Resolution',300);

pfinal = 0.5*( pnew(:,2:end) + pnew(:, 1:end-1));
pfinal = 0.5*( pfinal(2:end,:) + pfinal(1:end-1,:));

%%%%%%--------Pressure --------------------------------
figure
contourf(X,Y,pfinal,51,"LineStyle","none");
hold on
plot(XLnew(:,1), XLnew(:,2), "m",XRnew(:,1), XRnew(:,2), "m", "LineWidth",2.0);
plot(Xel(:,1), Xel(:,2), "k",Xer(:,1), Xer(:,2), "k", "LineWidth",2.0);
title("final pressure profile")
xlabel("eyewall ")
ylabel("p(x,y)")
colormap('jet')
xlim([-Lx,Lx]), ylim([0,Ly/5]);
colorbar
pp1= gcf;
exportgraphics(pp1,'Pressure.png','Resolution',300);

%%% --------Ploting the saccadic waveform ------------------
figure
subplot(2,2,1)
plot(time, Theta, 'b', 'LineWidth',1.2);
title('\bf Ang. displacement')
xlabel('\bf Time [sec]'); ylabel('\bf \theta(t) [deg]')
axis tight
grid on;

subplot(2,2,2)
plot(time, WallVel, 'b', 'LineWidth',1.2);
title('\bf Lin. Velocity')
xlabel('\bf Time [sec]'); ylabel('\bf v(t) [m/s]')
axis tight
grid on;

subplot(2,2,3)
plot(time, S, 'b', 'LineWidth',1.2);
title('\bf Lin. displacement')
xlabel('\bf Time [sec]'); ylabel('\bf S(t) [m]')
axis tight
grid on;

subplot(2,2,4)
plot(time, Omega, 'b', 'LineWidth',1.2);
title('\bf Ang. Velocity')
xlabel('\bf Time [sec]'); ylabel('\bf \omega(t) [deg/sec]')
axis tight
grid on;
pw= gcf;
exportgraphics(pw,'Saccade_waveform.png','Resolution',300);

%%%%%%%%%%%====================================================
%%%-------- Snapshots of the flow ---------------------------------
ksnap = floor((Num_SaveData + 1)/4 );
snapInd = ksnap: ksnap: Num_SaveData + 1;

figure
for js = 1:4
    kc = snapInd(js);
    subplot(2,2,js)
    contourf(X,Y,SaveMag(:,:,kc),"LineStyle","none");
    hold on
    plot(SaveXL(:,1,kc), SaveXL(:,2,kc), "r-",SaveXR(:,1,kc), SaveXR(:,2,kc), "r-","LineWidth",1.20);
    plot(Xel(:,1), Xel(:,2), "k-", Xer(:,1), Xer(:,2), "k-","LineWidth",0.80);
    sc = streamslice(X,Y, SaveU(:,:,kc), SaveV(:,:,kc), 6);
    set(sc,'LineWidth',0.05,'Color','w'); %[0.8, 0.55, 0.40]

    SimTime = sprintf('Time = %.2f s', SnapTime(kc) );
    text(0.5, 0.70, SimTime, 'Units','normalized','FontSize',10, 'color','w','FontWeight', ...
        'bold','BackgroundColor','k', 'HorizontalAlignment','center' );

    colorbar;

    grid on;
    xlim([-Lx,Lx]); ylim([0, Ly/2]);
    xlabel("\bf eyewall [m]"); ylabel("\bf vitreous [m]");
    colorbar
    % title('\bf Flow profile');
end
pfs= gcf;
exportgraphics(pfs,'Flow_Snapshot.png','Resolution',300);


%%%-------- Snapshot of the linker proteins ---------------------------------
figure
for jp = 1:4
    kp = snapInd(jp);
    subplot(2,2,jp)
    plot(Xel(:,1), SaveRho(:,1,1),'r', Xer(:,1), SaveRho(:,2,1), 'r','LineWidth',1.2);
    hold on
    plot(SaveXL(:,1,kp), SaveRho(:,1,kp),'b', SaveXR(:,1,kp), SaveRho(:,2,kp), 'b','LineWidth',1.2);

    SimTime = sprintf('Time = %.2f s', SnapTime(kp) );
    text(0.06*Lx, 0.80*rho_0init, SimTime, 'FontSize',10, 'color','w','FontWeight', ...
        'bold','BackgroundColor','k' );
   
    grid on;
    xlim([-Lx,Lx]);
    xlabel("\bf eyewall [m]"); ylabel("\bf \rho_b");
    % title('\bf Adhesion snapshot');
end
pas= gcf;
exportgraphics(pas,'Adhesion snapshot.png','Resolution',300);

%%%%%%%%%%==========Saving data into Excell ============================
rData = zeros(Nb+1, 4);
rData(2:end, :) = [XLnew(:,1), rhobL, XRnew(:,1), rhobR];    % Rhob_data
filename1 = 'RhobData.xlsx';                                 % Excel file name
writematrix(rData, filename1);                               % Writing results to Excel
headers = {'sL', 'Left_rhob','sR', 'Right_rhob' };           % Assigning headers
writecell(headers, filename1);


detachData = zeros(NumOfIter, 3);
detachData(2:end, :) = DLent;
filename2 = 'Detach_Length.xlsx';                      % Excel file name
writematrix(detachData, filename2);                    % Writing results to Excel
headers2 = {'Time', 'Left_Detach', 'Right_Detach' };   % Assigning headers
writecell(headers2, filename2);

end