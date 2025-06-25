function [Au,Av] = UV_LapMat(dt)

global Nx Ny h Re Uint Vint

%%% This function computes the Laplacian matrices Au & Av used to
%%% solve the Poisson equation involving u* & v* respectively.

Cx  = -(dt)/(2*Re*h^2);
Cxy = 1 + (2*dt/(Re*h^2) );
Cy  = -(dt)/(2*Re*h^2);

%%%%%%%==========================================================================================
%% ---Computing the laplacian matrix Au use to solve Poisson equation involving u*---------------
%%%%%%%==========================================================================================
upU = repmat([0 ; ones(Nx,1)], Ny, 1);
downU = repmat([ones(Nx,1) ; 0], Ny, 1);
MainU = repmat(ones(Nx+1,1), Ny, 1);

Au = spdiags([Cy*MainU, Cx*downU, Cxy*MainU, Cx*upU, Cy*MainU],[-(Nx+1), -1, 0, 1, Nx+1], Uint, Uint);
Au = full(Au);

Au(1:Nx+1, 1:Nx+1) = Au(1:Nx+1, 1:Nx+1) - Cy*eye(Nx+1);
Au(end-Nx:end, end-Nx:end) = Au(end-Nx:end, end-Nx:end) + Cy*eye(Nx+1);

for j = 1:Ny
    j1 = (j-1)*(Nx+1)+1;
    j2 = j*(Nx+1);
    Au(j1,j2) = Au(j1,j2)+ Cx;
    Au(j2,j1) = Au(j2,j1)+ Cx;
end


%%%%%%%==========================================================================================
%% ---Computing the laplacian matrix Av use to solve Poisson equation involving v*---------------
%%%%%%%==========================================================================================
upV = repmat([0 ; ones(Nx-1,1)], Ny-1, 1);
downV = repmat([ones(Nx-1,1) ; 0], Ny-1, 1);
MainV = repmat(ones(Nx,1), Ny-1, 1);
Av = spdiags([Cy*MainV, Cx*downV, Cxy*MainV, Cx*upV, Cy*MainV], [-(Nx), -1, 0, 1, Nx], Vint, Vint);
Av = full(Av);

for i = 1: Ny-1
    i1 = (i-1)*(Nx)+1;
    i2 = i*(Nx);
    Av(i1, i2) = Av(i1, i2) + Cx;
    Av(i2, i1) = Av(i2, i1) + Cx;
end

end