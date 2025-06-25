function [distL, distR] = DetachRate(wL, wR,rhoL0, rhoR0, Xl, Xl0, Xr, Xr0, dt)
%%% This function computes the length of detached retina at any time.
global rhoMin s

[rhobL,rhobR,~] = RhoB(wL, wR, rhoL0, rhoR0, Xl, Xl0, Xr, Xr0, dt);
%%%%------Left retina ------------------------------------
sl=  -s;                     %% x-component of the left retina
wLind =  rhobL >= rhoMin;    %% binary index for attached and detached retina {detached = 1, attached = 0}
lastIndxL = find(wLind == 1, 1, 'first');   %% finding the index of the clamped edge
distL = -sl(lastIndxL);      %% computing the length of detached retina

%%%%------Right retina ------------------------------------
sr = s;                  %% x-component of the right retina
wRind =  rhobR >= rhoMin;   %% binary index for attached and detached retina {detached = 1, attached =0}
lastIndxR = find(wRind == 1, 1, 'first');   %% finding the index of the clamped edge  
distR = sr(lastIndxR);    %% computing the length of detached retina


end