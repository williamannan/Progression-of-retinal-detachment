function [unew, vnew, XLnew, XRnew,pnew,rhobL, rhobR]=P1FirstIter(u0,v0,dt,Const_p, WallAcc,WallVel,Nt)
%%%% pp,pc, pu, pv
global h Nx Ny vS vN Xel Xer alpha

% time = 0:dt:tf;
%%%%%%---------User define functions ----------------------------------
[Ap] = PLapMat(Const_p);
[Xl, Xr] = LagStructure();         %% second Lagrangian structure
Xl0 = Xel;     Xr0 = Xer;
[rhoL0, rhoR0] = InitRhob(Xl, Xr);
[FaL0, FaR0] = Falpha0(u0, v0, Xl, Xl0, Xr, Xr0,dt);
[wL0, wR0] = BondDispW(Xl, Xr);

for step = 1: Nt+1
  
    ustar = zeros(Ny+2, Nx+1);
    vstar = zeros(Ny+1, Nx+2);
    pnew =  zeros(Ny+2, Nx+2);
  
    uBC = WallVel(2);

    [RHSu, RHSv] = P1RHS(FaL0, FaR0, u0,v0,Xl,Xl0,Xr,Xr0,dt);

    ustar(2:end-1, :) = RHSu;
    ustar(1,:) = 2*uBC - ustar(2, :);  
    ustar(end, :) = ustar(end-1,:);                                   %% Neumann at the top
    ustar(:,1) = ustar(:,end);     ustar(:,end) = ustar(:,1);         %% Periodic on the left and irght

    vstar(2:end-1, 2:end-1) = RHSv; 
    vstar(:, 1) = vstar(:, end-1);                           %% Periodic BC on the left
    vstar(:, end) = vstar(:, 2);                             %% Periodic BC on the left
    vstar(1,:) = vS;               vstar(end,:) = vN;          %% Dirichlet at top and bottom

    %%% Solving the Poisson equation for p
    dustar_dx = (ustar(2:end-1, 2:end) - ustar(2:end-1, 1:end-1) )/h;
    dvstar_dy = (vstar(2:end, 2:end-1) - vstar(1:end-1, 2:end-1) )/h;
    DivU = dustar_dx + dvstar_dy;

    vec_DivU = reshape(DivU', (Nx)*(Ny), 1 );
    p_int = Ap\((1/dt)*(vec_DivU));
    pnew(2:end-1, 2:end-1) = reshape(p_int, Nx, Ny)';

    pnew(1,:) = pnew(2,:);    pnew(end,:) = pnew(end-1, :);
    pnew(:,1) = pnew(:,2);    pnew(:,end) =  pnew(:,end-1);
    
    %%% Finding updated velocity
    dp_dx = (pnew(:,2:end) - pnew(:, 1:end-1))/h;
    dp_dy = (pnew(2:end, :) - pnew(1:end-1, :) )/h;

    unew = ustar - (dt) *dp_dx;
    unew(1,:) = 2*uBC -unew(2, :);
    unew(end, :) = unew(end-1,:);
    unew(:,1) = unew(:,end);      unew(:,end) = unew(:,1);

    vnew = vstar - (dt)*dp_dy;
    vnew(:, 1) = vnew(:, end-1);            %% Periodic BC on the left
    vnew(:, end) = vnew(:, 2); 
    vnew(1,:) = vS;     vnew(end,:) = vN;
    
   
    [XLnew, XRnew,~] = UpdateLag(WallAcc(2), wL0, wR0, rhoL0, rhoR0, FaL0, FaR0,u0,v0,Xl, Xl0, Xr,Xr0, dt);
    

    %% ------ Update initial data ----------------
    [UibL,UibR] = Uib(unew,vnew, XLnew, XRnew);

    [rhobL,rhobR,~] = RhoB(wL0, wR0, rhoL0, rhoR0, Xl, Xl0, Xr, Xr0, dt);
    
    dxLdt = (XLnew - Xl)/dt;
    dxRdt = (XRnew - Xr)/dt;

    FaL0 = FaL0 + alpha*(UibL - dxLdt )*dt;
    FaR0 = FaR0 + alpha*(UibR - dxRdt )*dt;
    [wL, wR] = BondDispW(XLnew, XRnew);

    rhoL0 = rhobL; 
    rhoR0 = rhobR;

    u0 = unew;
    v0 = vnew;

    Xl0 = Xl;     Xl = XLnew;

    Xr0 = Xr;     Xr = XRnew;

    wL0 = wL;     wR0 = wR;

end

end