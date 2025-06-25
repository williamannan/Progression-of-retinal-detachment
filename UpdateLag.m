function [XLnew, XRnew, cfl]=UpdateLag(WallAcc, wL, wR, rhoL0, rhoR0, FaL0, FaR0, u,v,Xl, Xl0, Xr,Xr0, dt)

global ds Nb Xel Xer Kend

[TL, TR] = Tension(WallAcc, wL, wR, rhoL0, rhoR0, FaL0, FaR0, u, v, Xl,Xl0, Xr,Xr0,dt);  %% Tension
[FhL, FhR,cfl] = AdhesionForce(wL, wR, rhoL0, rhoR0, Xl, Xl0, Xr, Xr0, dt);
[FaL, FaR] = LagForce(FaL0, FaR0, u, v, Xl, Xl0, Xr, Xr0,dt);           %% Force on the material
[FbL, FbR] = BendForce(Xl, Xr);
[AtL]=TentionMat(TL);
[AtR]=TentionMat(TR);
I = eye(Nb);

%%%%%%%%%=======================================================
ALmat = I - (dt^2)*AtL;      

glvec = zeros(Nb, 2);
glvec(1,:) = 2*Kend*norm((Xr(1,:) - Xl(1,:)), 2)*( Xl(2, :) - Xl(1,:) )/(ds^2);
glvec(end,:) = WallAcc + FaL(end, :) + FhL(end, :) - FbL(end, :);  %% Term coming from Fl

rhsL = (2*Xl - Xl0) + (dt^2)*(FbL - FaL - FhL + glvec);
XLnew = ALmat\rhsL;
XLnew(end, :) = Xel(end, :);

%%%%%%%%%-------------------------------------------
ARmat = I - (dt^2)*AtR;
grvec = zeros(Nb, 2);
grvec(1,:) =   2*Kend*norm((Xr(1,:) - Xl(1,:)), 2)*( Xr(2, :) - Xr(1,:) )/(ds^2); 
grvec(end, :)  = WallAcc + FaR(end, :) + FhR(end, :) -FbR(end, :) ;   %% Term coming from Fl

rhsR = (2*Xr - Xr0) + (dt^2)*(FbR - FaR- FhR + grvec);
XRnew = ARmat\rhsR;
XRnew(end,:) = Xer(end, :);


%%%%%%%%%%%%%%%%===============================================
      function [ At]=TentionMat(T)
        mdiag = [T; 0];
        mdiag(2:end-1) = mdiag(2:end-1) + T(1:end-1);
        updiag = [0; T];
        lowdiag = [T; 0];

        At = spdiags([lowdiag, -mdiag, updiag], [-1, 0, 1], Nb, Nb );
        % At = full(At);
        At(1, :) = 2*At(1, :);
        At(end, :) = 0;
        At = (1/(ds^2))*At;
    end

end