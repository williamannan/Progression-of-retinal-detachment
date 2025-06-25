function [dh] = Delta(r)
global h
if (r >= -2 && r <= -1)
    dh = 1/(8*h)*( 5 + 2*r - sqrt(-7 - 12*r - 4*r^2) );
elseif (r >= -1 && r <= 0) 
    dh = 1/(8*h)*( 3 + 2*r + sqrt(1 - 4*r - 4*r^2) );
elseif (r >= 0 && r <= 1)
    dh = 1/(8*h)*( 3 - 2*r + sqrt(1 + 4*r - 4*r^2) );
elseif (r>=1 && r<= 2)
    dh = 1/(8*h)*(5 - 2*r - sqrt(-7 + 12*r - 4*r^2) );
else
    dh = 0;
end

end