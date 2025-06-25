function [RHSu, RHSv]= P2_RHS(FaL0, FaR0,u0,un,v0,vn,p0,Xl,Xl0, Xr,Xr0,dt, uBC)

global Nx Ny h Re vS vN Usize Vsize

% [RHSu,RHSv] = Grad2U(u0,v0);
[uDu0,vDv0] = AdvecTerm(u0,v0);
[uDun,vDvn] = AdvecTerm(un,vn);

Hu = (3/2)*uDun - (1/2)*(uDu0);
Hv = (3/2)*vDvn - (1/2)*(vDv0);

[LapU,LapV] = Grad2U(un,vn);

[Px,Py] = GradP(p0);

[Bu,Bv] = BC_Func(dt, uBC);

[Fu, Fv] = FluidForce(FaL0, FaR0, un,vn, Xl,Xl0, Xr,Xr0, dt);

RHSu = un(2:end-1, :) - dt*Hu- dt*Px + (dt)/(2*Re)*LapU + dt*Fu(2:end-1,:) + Bu;
RHSv = vn(2:end-1, 2:end-1)- dt*Hv - dt*Py+ (dt)/(2*Re)*LapV + dt*Fv(2:end-1, 2:end-1) + Bv;

%%% Converting the RHS into  vector;
% Usize = (Ny)*(Nx+1);
% Vsize = (Ny-1)*(Nx);
RHSu = reshape(RHSu', Usize, 1);

RHSv = reshape(RHSv', Vsize, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%============================================================================
    function [uDu,vDv] = AdvecTerm(u,v)
        %%%%%%%%%*********************************************************************************
        %% ------------ computing the advection term uDu = u*u_x + v*u_y ----------------------------
        %%% ------------------- computing u*u_x -------------------------------
        UxRight = [1:Nx+1, 1];
        UxLeft = [Nx+1, 1:Nx+1];
        uh_avg = 0.5*( u(:, UxRight) + u(:, UxLeft) );
        ux = ( u(:, UxRight) - u(:, UxLeft) )/h;

        uux_prod = ux.*uh_avg;
        uux = 0.5*(uux_prod(:,2:end) + uux_prod(:,1:end-1) );    %% Averaging the product to get u*u_x

        %%% ------------------- computing v*u_y -------------------------------
        vh_avg  = 0.5*( v(:,2:end) + v(:,1:end-1) );             %% horizontal averaging of v
        uvet_diff = ( u(2:end,:) - u(1:end-1,:) )/h;            %% certical differencing of u
        vuy_prod = vh_avg.*uvet_diff;                            %% product oof the above two
        vuy = 0.5*( vuy_prod(2:end,:) + vuy_prod(1:end-1,:) );   %% Averaging the product to get v*u_y

        %%% --------- computing the final output uDu = u*u_x + v*u_y--------------------
        uDu = uux(2:end-1,:) + vuy;

        %%%%%%%%%*********************************************************************************
        %% ------------------ computing the advection term vDv = u*v_x + v*v_y -----------------------
        %%% -------------------- Computing u*v_x -----------------------------------
        uvert_avg = 0.5*( u(2:end,:) + u(1:end-1,:) );           %% vertical averaging of u
        vh_diff = ( v(:,2:end) - v(:, 1:end-1) )/h;              %% horizontal difference of v
        uvx_prod = uvert_avg.*vh_diff;                           %% product oof the above two
        uvx = 0.5*(uvx_prod(:,2:end) + uvx_prod(:,1:end-1));     %% Averaging the product to get u*v_x

        %%% -------------------- Computing v*v_y -----------------------------------
        vvert_avg = 0.5*( v(2:end,:) + v(1:end-1,:) );           %% vertical averaging of v
        vvert_diff = ( v(2:end,:) - v(1:end-1,:) )/h;            %% vertical differencing of v
        vvy_prod = vvert_avg.*vvert_diff;                        %% product oof the above two
        vvy = 0.5*( vvy_prod(2:end,:) + vvy_prod(1:end-1,:));    %% Averaging the product to get v*v_y

        %%% computing the final output vDv = u*v_x + v*v_y
        vDv = ( uvx(2:end-1,:) + vvy(:,2:end-1));

    end     %% end of advection term function

%%%%%%%==================================================================================
 %%%%%%%==================================================================================
    function [LapU, LapV] = Grad2U(u,v)
        %%% Compute the laplacian grad^2 u = u_xx + u_yy and grad^2 v = v_xx + v_yy
      
        UxxRight = [2:Nx+1, 1];
        UxxLeft = [Nx+1, 1:Nx];
        %% -------------- LapU = u_xx + u_yy -----------------------------------------------------
       
        u_xx = ( u(:,UxxRight) - 2*u + u(:, UxxLeft) )/(h^2);


        u_yy = ( u(3:end, :) -2*u(2:end-1, :) + u(1:end-2, :))/(h^2);

        LapU = u_xx(2:end-1, :) + u_yy;

        %%%%%%%%%*********************************************************************************
        %% -------------- LapV = v_xx + v_yy -----------------------------------------------------
        v_xx = ( v(:, 3:end) - 2*v(:, 2:end-1) + v(:,1:end-2) )/(h^2);
        v_yy = ( v(3:end, :) -2*v(2:end-1, :) + v(1:end-2, :) )/(h^2);
        LapV = v_xx(2:end-1, :) + v_yy(:,2:end-1);
    end

%%%%%%%==================================================================================
    function [Px,Py] = GradP(p)
        %% ------------ Computing Gradiend p --------------------
        dp_dx = ( p(:, 2:end) - p(:, 1:end-1))/h;
        Px = dp_dx(2:end-1, :);

        dp_dy = (p(2:end , :) - p(1:end-1, :) )/h;
        Py = dp_dy(2:end-1 , 2:end-1);
    end   %% end of gradient function fo

%%%%%%%==========================================================================================
    function [Bu,Bv] = BC_Func(dt, uBC)
        %% -Organizing the boundary terms from the discretization of  the diffusion term moved to the RHS
        %%% Refer to discretization details and terms moved to the RHS for what is done here.
        Cy = -(dt)/(2*Re*h^2);         %%Coefficient from discretization of diffusion eqn

        Bu = zeros(Ny, Nx+1); 
        Bu(1,:) = -2*Cy*uBC;

        Bv = zeros(Ny-1, Nx);             %% Number of interior nodes for
        Bv(1,:) = -Cy*vS;
        Bv(end,:) = -Cy*vN;

    end         %% end of function organizing the BC moved to the right

end