function [rhoL0, rhoR0] = InitRhob(Xl, Xr)

global rho0 eta lx0

rhoL0 = rho0*( 1 + exp(eta*(Xl(:,1) + lx0)) ).^(-1);  
rhoR0 = rho0*( 1 + exp(-eta*(Xr(:,1) - lx0)) ).^(-1); 


end