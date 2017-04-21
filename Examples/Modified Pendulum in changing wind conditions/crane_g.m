function gx = crane_g( x,t,var )
%% CRANE_G

% Retrieve constants
m = var(1);
g = var(2);

alpha = 1 - 0.04*t;
l = 0.5 + 0.03*(50-t);

% Evaluate function
gx = zeros(3,1);

gx(1) = x(1)^2 + x(2)^2 - l^2;
gx(2) = 2*x(1)*x(3) + 2*x(2)*x(4);
gx(3) = 2*x(3)^2 + 2*x(4)^2 + 2*alpha*x(1) - 4/m*(x(1)^2 + x(2)^2)*x(5) - 2*g*x(2);