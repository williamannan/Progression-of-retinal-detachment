function [TL, TR] = Tension(WallAcc, wL, wR, rhoL0, rhoR0, FaL0, FaR0, u, v, Xl,Xl0, Xr,Xr0,dt)

global ds Nb Kend

[FbL, FbR] = BendForce(Xl, Xr);
[FaL, FaR] = LagForce(FaL0, FaR0, u, v, Xl, Xl0, Xr, Xr0,dt);
[FhL, FhR, ~] = AdhesionForce(wL, wR, rhoL0, rhoR0, Xl, Xl0, Xr, Xr0, dt);
AtL = Tmat(Xl);                %% LHS tension matrix
AtR = Tmat(Xr);                %% RHS tension matrix

%%%%%%%%-----------------------------------------------------------------------------------
Yl = ( Xl(2:end, :) - Xl(1:end-1, :) )/ds;            %%% DsX_n  for the left retina
Yl0 = ( Xl0(2:end, :) - Xl0(1:end-1, :) )/ds;         %%% DsX_n-1 for the left retina

Uln = (Xl - Xl0)/dt;                                    %% Velocity D_t(X) = (X_n - X_n-1)/dt 
DsUln = (Uln(2:end, :) - Uln(1:end-1, :))/ds;           %% D_s(u) = (u_j+1 - u_j)/ds

%%% --------------------- RHS Left retina--------------------------
RHSl1 = ( 1 - 2*sum(Yl.*Yl, 2) +  sum(Yl0.*Yl0, 2) )/(2*dt^2);

RHSl2 = sum(DsUln.*DsUln, 2);

DsXl = (Xl(2:end, :) - Xl(1:end-1, :))/ds;              %% Ds(X*) = (X*_j+1 - X*_j)/ds
Flterm = FbL - FaL - FhL;                               %% Fb - F
DsFlterm = (Flterm(2:end, :) - Flterm(1:end-1, :))/ds;    %% Ds(Fb - F)

RHSl3 = sum(DsXl.*DsFlterm, 2);                               %% (DsX*)*[Ds(Fb - F)]

LHSl1 = zeros(Nb-1, 1);
LHSl1(1) = -2*Kend*norm((Xr(1, :) - Xl(1, :)), 2)*sum(Yl(1, :).*Yl(1, :))/(ds^2);
LHSl1(end) =sum(DsXl(end, :).*( (FaL(end, :) + FhL(end, :) - FbL(end, :) + WallAcc))/ds );
TL = AtL\( RHSl1 - RHSl2 - RHSl3 - LHSl1);


%%%%%%%%-----------------------------------------------------------------------------------
Yr = ( Xr(2:end, :) - Xr(1:end-1, :) )/ds;            %%% DsX_n  for the left retina
Yr0 = ( Xr0(2:end, :) - Xr0(1:end-1, :) )/ds;         %%% DsX_n-1 for the left retina

Urn = (Xr - Xr0)/dt;                                  %% Velocity D_t(X) = (X_n - X_n-1)/dt 
DsUrn = (Urn(2:end, :) - Urn(1:end-1, :))/ds;         %% D_s(u) = (u_j+1 - u_j)/ds

%%% --------------------- RHS Left retina--------------------------
RHSr1 = ( 1 - 2*sum(Yr.*Yr, 2) +  sum(Yr0.*Yr0, 2) )/(2*dt^2);

RHSr2 = sum(DsUrn.*DsUrn, 2);

DsXr = (Xr(2:end, :) - Xr(1:end-1, :))/ds;     %% Ds(X*) = (X*_j+1 - X*_j)/ds
Frterm = FbR - FaR - FhR;                                       %% Fb - F
DsFrterm = (Frterm(2:end, :) - Frterm(1:end-1, :))/ds;    %% Ds(Fb - F)

RHSr3 = sum(DsXr.*DsFrterm, 2);                               %% (DsX*)*[Ds(Fb - F)]

LHSr1 = zeros(Nb-1, 1);
LHSr1(1) = -2*Kend*norm((Xr(1, :) - Xl(1, :)), 2)*sum(Yr(1, :).*Yr(1, :))/(ds^2);
LHSr1(end) = sum( DsXr(end, :).*( (FaR(end, :) + FhR(end, :) - FbR(end, :) + WallAcc))/ds );
TR = AtR\( RHSr1 - RHSr2 - RHSr3 - LHSr1);

%%%%%%%%%%%%%%%%%%%%%===============================================================
    function [At] = Tmat(X)
        Y = (X(2:end, :) - X(1:end-1, :) )/ds;
        mdiag = sum(Y.*Y, 2);    
        updiag = [0 ; sum(Y(1:end-1, :).*Y(2:end, :), 2)];
        lowdiag = [sum(Y(1:end-1, :).*Y(2:end, :), 2) ; 0];

        At = spdiags([lowdiag, -2*mdiag, updiag], [-1, 0, 1], Nb-1, Nb-1 );
        At(1,1) =  (3/2)*At(1,1);
        At(end, end) = (1/2)*At(end, end); 
        At = (1/ds^2)*At;
    end

end