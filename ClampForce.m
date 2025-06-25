function [CFtot,Ffl] = ClampForce(WallAcc, wL, wR, rhobL,rhobR,rhoL0, rhoR0, FaL0, FaR0, u, v, Xl,Xl0, Xr,Xr0,dt)
global ds Nb Kend rhoMin

[TL, TR] = Tension(WallAcc, wL, wR, rhoL0, rhoR0, FaL0, FaR0, u, v, Xl,Xl0, Xr,Xr0,dt);
[FbL, FbR] = BendForce(Xl, Xr);
[FaL, FaR] = LagForce(FaL0, FaR0, u, v, Xl, Xl0, Xr, Xr0,dt);
[FhL, FhR, ~] = AdhesionForce(wL, wR, rhoL0, rhoR0, Xl, Xl0, Xr, Xr0, dt);
%%% =====  Computing the tension force =========================
%%% ----- Left retina --------------
DsXL = (Xl(2:end, :) - Xl(1:end-1, :))/ds;                          %% dX/ds
TDsXL = TL.*DsXL;                                                   %% T*dX/ds
FtL = zeros(Nb, 2);                                                  
FtL(2:end-1, :) = (TDsXL(2:end, :) - TDsXL(1:end-1, :))/ds;
TDsXL0 = Kend*norm((Xr(1,:) - Xl(1,:)), 2)*(Xl(2,:) - Xl(1, :))/ds;
FtL(1,:) = (TL(1)*DsXL(1,:) - TDsXL0 )/(ds/2);
FtL(end, :) = WallAcc - FbL(end, :) + FaL(end, :) + FhL(end, :);

TotalFL = FtL + FbL - FaL - FhL;
FLdrive = FaL;

LeftClamInd =  rhobL >= rhoMin;  
lastIndxL = find(LeftClamInd == 1, 1, 'first');

%%% ----- right retina --------------
DsXR = (Xr(2:end, :) - Xr(1:end-1, :))/ds;
TDsXR = TR.*DsXR;
FtR = zeros(Nb, 2);
FtR(2:end-1, :) = (TDsXR(2:end, :) - TDsXR(1:end-1, :))/ds;
TDsXR0 = Kend*norm((Xr(1,:) - Xl(1,:)), 2)*(Xr(2,:) - Xr(1, :))/ds;
FtR(1,:) = (TR(1)*DsXR(1,:) - TDsXR0 )/(ds/2);
FtR(end, :) = WallAcc - FbR(end, :) + FaR(end, :) + FhR(end, :);

TotalFR = FtR + FbR - FaR -FhR;
FRdrive = FaR; 

RightClamInd =  rhobR >= rhoMin;  
lastIndxR = find(RightClamInd == 1, 1, 'first');

% ClampFtot = [TotalFL(lastIndxL,2),TotalFR(lastIndxR,2)];
% ClampFeffect = [FLdrive(lastIndxL,2), FRdrive(lastIndxR,2)];

CFtot = [abs(TotalFL(lastIndxL,2)),abs(TotalFR(lastIndxR,2))];
Ffl = [abs(FLdrive(lastIndxL,2)), abs(FRdrive(lastIndxR,2))];

end