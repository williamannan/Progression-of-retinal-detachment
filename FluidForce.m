function [Fu, Fv] = FluidForce(FaL0, FaR0, u,v, Xl,Xl0, Xr,Xr0, dt)

global Nx Ny h Nb ds xu xv yu yv rho
%%% Note: Always, (Xl, Xr) represent current value while  (Xl0, Xr0) represent previous value. 

[FaL, FaR] = LagForce(FaL0, FaR0, u, v, Xl, Xl0, Xr, Xr0,dt);

%% -------- Computing the fluid force --------------
%%% u component of the force exert on the fluid by the immersed boundary
flu = zeros(Ny+2, Nx+1);         %% Initializing the u component of the fluid force on the left retina
fru = zeros(Ny+2, Nx+1);         %% Initializing the u component of the fluid force on the right retina

for iy = 1: Ny+2
    for ix = 1:Nx+1
        smlu = 0;
        smru = 0;
        for rs = 1:Nb
            l1 = (xu(ix) - Xl(rs,1))/h;        %% (x - X)
            l2 = (yu(iy) - Xl(rs, 2))/h;       %% (y - Y)
            smlu = smlu + FaL(rs,1)*Delta(l1)*Delta(l2)*ds;

            r1 = (xu(ix) - Xr(rs,1))/h;        %% (x - X)
            r2 = (yu(iy) - Xr(rs, 2))/h;       %% (y - Y)
            smru = smru + FaR(rs,1)*Delta(r1)*Delta(r2)*ds;
        end
        flu(iy, ix) = smlu;
        fru(iy, ix) = smru;
    end
end
Fu = rho*(flu + fru);       %% total u component of force on the fluid

%%% v component of the force exert on the fluid by the immersed boundary
flv = zeros(Ny+1, Nx+2);            %% Initializing the v component of the fluid force on the left retina
frv = zeros(Ny+1, Nx+2);            %% Initializing the v component of the fluid force on the right retina

for jy = 1:Ny+1
    for jx = 1: Nx+2
        smlv = 0;
        smrv = 0;
        for ks = 1:Nb
            %%% force on the right half of detached retina
            kl1 = (xv(jx) - Xl(ks, 1))/h;
            kl2 = (yv(jy) - Xl(ks, 2))/h;
            smlv = smlv + FaL(ks, 2)*Delta(kl1)*Delta(kl2)*ds;

            kr1 = (xv(jx) - Xr(ks, 1))/h;
            kr2 = (yv(jy) - Xr(ks, 2))/h;
            smrv = smrv + FaR(ks, 2)*Delta(kr1)*Delta(kr2)*ds;
        end
        flv(jy, jx) = smlv;
        frv(jy, jx) = smrv;
    end
end
Fv = rho*(flv + frv);               %% total v component of force on the fluid

%%%%%------ If you want to visualize the force uncomment this part --------
% figure
% imagesc(xu, xv, Fu);
% title('u-compoenent force on fluid')
% xlabel("x"); ylabel('y');
% colorbar;
% 
% figure
% imagesc(xv, xv, Fv);
% title('v-compoenent force on fluid')
% xlabel("x"); ylabel('y');
% colorbar;

end