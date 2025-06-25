function [FaL0, FaR0] = Falpha0(u, v, Xl, Xl0, Xr, Xr0,dt)
global alpha 

%%% Note: Always, (Xl, Xr) represent current value while  (Xl0, Xr0) represent previous value. 

[UibL,UibR] = Uib(u,v, Xl, Xr);

dxLdt = (Xl - Xl0)/dt;
dxRdt = (Xr - Xr0)/dt;

FaL0 = alpha*(UibL - dxLdt )*dt;
FaR0 = alpha*(UibR - dxRdt )*dt;

end