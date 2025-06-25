function [FbL, FbR] = BendForce(Xl, Xr)
global gamma ds Nb

%%%%%%%%%%%%%%========= Using Extrapolation ================================
XLext = zeros(Nb+3, 2);                          %% Extended left retina for computing the 4th derivative
XLext(1,:) = 2*Xl(1, :) - Xl(2, :);              %% Imposing the boundary condition DssX|_0 = 0
XLext(2:end-2,:) = Xl;
XLext(end-1,:) = Xl(end,:) + ds*[-1,0];          %% Imposing the BC DsX = (-1, 0) for the left retina .

sl0_end = XLext(end-5:end-1, 1)';                 %% Last 5 Lagragian position  

sl_end = sl0_end(end) - ds;                       %% Last  ghost Lagrangian position

XLext(end, :) = InterpPoly(sl0_end, XLext(end-5:end-1, :), sl_end); %% Interpolating the last ghost point

D4XL = zeros(Nb, 2);
D4XL(1, :) = 2*( Xl(3,:) - 2*Xl(2,:) + Xl(1, :) )/(ds^4);
D4XL(2:end, :) = (XLext(5:end, :) - 4*XLext(4:end-1, :) + 6*XLext(3:end-2, :) ...
               - 4*XLext(2:end-3, :) + XLext(1:end-4, :) )/(ds^4);
FbL = -gamma*D4XL;

%%%%%%%%%%%%*******************---Right Retina -------******************************************
XRext = zeros(Nb+3, 2);                  %% Extended left retina for computing the 4th derivative
XRext(1,:) = 2*Xr(1, :) - Xr(2, :);      %% Imposing the boundary condition DssX|_0 = 0
XRext(2:end-2,:) = Xr;
XRext(end-1,:) = Xr(end,:) + ds*[1,0];
  
sr0_end = XRext(end-5:end-1, 1)';        %% Last 5 Lagragian position  

sr_end = sr0_end(end) + ds;              %% Last 2 ghost Lagrangian position

XRext(end, :) = InterpPoly(sr0_end, XRext(end-5:end-1, :), sr_end);

D4XR = zeros(Nb, 2);
D4XR(1, :) = 2*(Xr(3,:) - 2*Xr(2,:) + Xr(1, :))/(ds^4);
D4XR(2:end, :) = (XRext(5:end, :) - 4*XRext(4:end-1, :) + 6*XRext(3:end-2, :) ...
               - 4*XRext(2:end-3, :) + XRext(1:end-4, :) )/(ds^4);
FbR = -gamma*D4XR;

end