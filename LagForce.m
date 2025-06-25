function [FaL, FaR] = LagForce(FaL0, FaR0, u, v, Xl, Xl0, Xr, Xr0,dt)

global beta Kend

%%% Note: Always, (Xl, Xr) represent current value while  (Xl0, Xr0) represent previous value.
[UibL,UibR] = Uib(u,v, Xl, Xr);

dxdt_L = ( Xl - Xl0)/dt;
dxdt_R = ( Xr - Xr0)/dt;


FaL = FaL0 + beta*(UibL  - dxdt_L);
FaR = FaR0 + beta*(UibR  - dxdt_R);

%%% ---------- Force at the free end -----------------
FaL(1,:) = FaL(1,:) - Kend*(Xr(1,:) - Xl(1,:));
FaR(1,:) = FaR(1,:) + Kend*(Xr(1,:) - Xl(1,:));

end