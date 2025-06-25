function [wL, wR] = BondDispW(Xl, Xr)

global w0
wL = (Xl(:,2) - w0);       % vertical displacement of the left retina;
wR = (Xr(:,2) - w0);       % vertical displacement of the right retina;
end