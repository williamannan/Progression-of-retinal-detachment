function [RHSu, RHSv] = P1RHS(FaL0, FaR0, u,v,Xl,Xl0, Xr,Xr0, dt)

global Nx h Re

[Fu, Fv] = FluidForce(FaL0, FaR0, u,v, Xl,Xl0, Xr,Xr0, dt);

[uDu,vDv] = AdvecTerm(u,v);

[LapU, LapV] = Grad2U(u,v);

RHSu = u(2:end-1, :) - dt*uDu + (dt/Re)*LapU + dt*Fu(2:end-1,:);
RHSv = v(2:end-1, 2:end-1) - dt*vDv + (dt/Re)*LapV + dt*Fv(2:end-1,2:end-1);

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

end