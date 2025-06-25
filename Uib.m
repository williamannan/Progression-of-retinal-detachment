function [UibL,UibR] = Uib(u,v, Xl, Xr)
global Nx Ny h Nb xu xv yu yv

%%% ---  X-component of Lag position --------
uXl = zeros(Nb,1);         %% Initializing the x-component of the left immersed boundary
uXr = zeros(Nb,1);         %% Initializing the x-component of the right immersed boundary

for k = 1: Nb
    smlu = 0;
    smru = 0;
    for iy = 1:Ny+2
        for ix = 1:Nx+1
            l1 = (xu(ix) - Xl(k, 1))/h;
            l2 = (yu(iy) - Xl(k, 2) )/h;
            smlu = smlu + u(iy, ix)*Delta(l1)*Delta(l2)*h*h;

            r1 = (xu(ix) - Xr(k, 1))/h;
            r2 = (yu(iy) - Xr(k, 2) )/h;
            smru = smru + u(iy, ix)*Delta(r1)*Delta(r2)*h*h;
        end
    end
    uXl(k) = smlu;
    uXr(k) = smru;
end

%%% ---  Y-component of Lag position --------
uYl = zeros(Nb,1);     %% Initializing the y-component of the left immersed boundary
uYr = zeros(Nb,1);     %% Initializing the y-component of the right immersed boundary

for kk = 1: Nb
    smlv = 0;
    smrv = 0;
    for jy = 1:Ny+1
        for jx = 1:Nx+2
            ll1 = (Xl(kk, 1) - xv(jx) )/h;  %(xv(jx) - Xl(kk, 1))/h;
            ll2 = (Xl(kk, 2) - yv(jy) )/h;  %(yv(jy) - Xl(kk, 2) )/h;
            smlv = smlv + v(jy, jx)*Delta(ll1)*Delta(ll2)*h*h;

            rr1 = (Xr(kk, 1) - xv(jx) )/h;  % (xv(jx) - Xr(kk, 1))/h;
            rr2 = (Xr(kk, 2) - yv(jy) )/h;  % (yv(jy) - Xr(kk, 2) )/h;
            smrv = smrv + v(jy, jx)*Delta(rr1)*Delta(rr2)*h*h;
        end
    end
    uYl(kk) = smlv;
    uYr(kk) = smrv;
end

UibL = [uXl, uYl];
UibR = [uXr, uYr];

end