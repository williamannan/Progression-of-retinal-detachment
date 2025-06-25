function [FhL, FhR, cfl] = AdhesionForce(wL, wR, rhoL0, rhoR0, Xl, Xl0, Xr, Xr0, dt)
global Cb Nb CbFactor s w0   % Xel Xer

XLrest = [-s', w0*ones(Nb, 1)];
XRrest = [s', w0*ones(Nb, 1)];

[rhobL,rhobR, cfl] = RhoB(wL, wR, rhoL0, rhoR0, Xl, Xl0, Xr, Xr0, dt);
nu = zeros(Nb, 2);               %% Initializing the the unit vector in the vertical direction
nu(:,2) = 1;                     %% Unit vector in the vertical direction

%%%------------------Left retina ------------------------------
ClPlus = (Xl(:,2) - XLrest(:,2)) >=0;
ClMinus = (Xl(:,2) - XLrest(:,2)) < 0;

CbL = Cb*(ClPlus + CbFactor*ClMinus);
FhL = (CbL.*rhobL).*(nu.*(Xl-XLrest));

%%%------------------Right retina ------------------------------
CrPlus = (Xr(:,2) - XRrest(:,2)) >=0;
CrMinus = (Xr(:,2) - XRrest(:,2)) < 0;

CbR = Cb*(CrPlus + CbFactor*CrMinus);

FhR = (CbR.*rhobR).*(nu.*(Xr-XRrest));

end