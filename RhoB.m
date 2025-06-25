function [rhobL, rhobR,cfl] = RhoB(wL, wR, rhoL0, rhoR0, Xl, Xl0, Xr, Xr0, dt)

global Xib Xiu sigmau sigmab rho0 Nb ds

wL = abs(wL);               %% Absolute value of the  left displacement
wR = abs(wR);               %% Absolute value of the right displacement

[DsXl, DsXr] = dX_ds(Xl, Xr);

%%%%%%%%%%%%%%%%%%%====================================================================
%% -------------- Left Retina --------------------------------------------------------
tau2L = sum(DsXl.*DsXl, 2);    %% DsX.DsX
Uln = (Xl - Xl0)/dt;
XLbar = sum(DsXl.*Uln, 2)./tau2L;

C1L = -(dt/ds)*XLbar;
C2L = 1 + (dt/ds)*XLbar + dt*( Xib*exp(-sigmab*wL) + Xiu*exp(sigmau*wL) );
ALmat = RrhobMat(C1L, C2L);
full(ALmat);
rhsL = rhoL0 + dt*rho0*Xib*exp(-sigmab*wL);
rhobL = ALmat\rhsL;
rhobL(1) = rhoL0(1);

%%%%%%%%%%%%%%%%%%%====================================================================
%% -------------- Right Retina --------------------------------------------------------
tau2R = sum(DsXr.*DsXr, 2);    %% DsX.DsX
Urn = (Xr - Xr0)/dt;
XRbar = sum(DsXr.*Urn, 2)./tau2R;

C1R = -(dt/ds)*XRbar;
C2R = 1 + (dt/ds)*XRbar + dt*( Xib*exp(-sigmab*wR) + Xiu*exp(sigmau*wR) );
ARmat = RrhobMat(C1R, C2R);
full(ARmat);
rhsR = rhoR0 + dt*rho0*Xib*exp(-sigmab*wR);
rhobR = ARmat\rhsR;
rhobR(1) = rhoR0(1);

cfl = (dt/ds)*[max(XLbar), max(XRbar)];

%%%%%%%%%%%%%%%%%%%====================================================================
%% -------------- User defined functions --------------------------------------------------------
    function [DsXl, DsXr] = dX_ds(XL, XR)
        DsXl = zeros(Nb, 2);
        DsXl(1, :) = (XL(2, :) - XL(1, :))/ds;
        DsXl(2:end-1, :) = (XL(3:end, :) - XL(1:end-2, :))/(2*ds);
        DsXl(end, :) = [-1, 0];

        DsXr = zeros(Nb, 2);
        DsXr(1, :) = (XR(2, :) - XR(1, :))/ds;
        DsXr(2:end-1, :) = (XR(3:end, :) - XR(1:end-2, :))/(2*ds);
        DsXr(end, :) = [1, 0];
    end


    function [Rmat] = RrhobMat(C1, C2)
        Rmat = spdiags([C1, C2], [-1, 0], Nb, Nb);
        Rmat(1,:) = 0;
        Rmat(1,1) = 1;
    end




end
