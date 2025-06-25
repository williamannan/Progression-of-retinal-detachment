function [Ap] = PLapMat(Const_p)

global Nx Ny h
Psize = (Ny)*(Nx);       %% size of the interior pressure nodes

Dx = 1/(h^2);
Dxy = -4/(h^2);
Dy = 1/(h^2);

upP = repmat([0 ; ones(Nx-1,1)], Ny, 1);
downP = repmat([ones(Nx-1,1) ; 0], Ny, 1);
MainP = repmat(ones(Nx,1), Ny, 1);
Ap = spdiags([Dy*MainP, Dx*downP, Dxy*MainP, Dx*upP, Dy*MainP], [-(Nx),-1, 0, 1, Nx], Psize, Psize);
Ap = full(Ap);
Ap(1:Nx,1:Nx) = Ap(1:Nx,1:Nx) + eye(Nx)*Dy;
Ap(end - Nx+1:end, end-Nx+1:end) = Ap(end - Nx+1:end, end-Nx+1:end) + eye(Nx)*Dy;

for i = 1:Ny
    j1 = (i-1)*(Nx) + 1;
    j2 = i*(Nx);
    Ap(j1, j1) = Ap(j1,j1) + Dx;
    Ap(j2, j2) = Ap(j2,j2) + Dx;
end
Ap(1,:) = 0;
Ap(1,1) = Const_p;


end