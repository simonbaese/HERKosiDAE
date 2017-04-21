function gx = massspringchain_g( x,t,var )
%% MASSSPRINGCHAIN_G

% Retrieve constants
m = var(1);
C = var(2);

% Evaluate function
gx = zeros(5,1);

gx(1) = x(2) - sin(t);
gx(2) = x(5) - cos(t);
gx(3) = C/m*(x(1) - x(2)) - C/m*(x(2) - x(3)) + sin(t);
gx(4) = C/m*(x(4) - x(5)) - C/m*(x(5) - x(6)) + cos(t);
gx(5) = C*(-3*C*(x(1)-x(2))+3*C*(x(2)-x(3))+2*x(7))/m^2 - sin(t);