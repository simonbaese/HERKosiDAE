function gx = pendulum_g( x,~,var )
%% PENDULUM_G

% Retrieve constants
m = var(1);
l = var(2);
g = var(3);

% Evaluate function
% gx = zeros(1,1);
% gx(1) = x(3)^2 + x(4)^2 - 2/m*(x(1)^2 + x(2)^2)*x(5) - g*x(2);

gx = zeros(3,1);
gx(1) = x(1)^2 + x(2)^2 - l^2;
gx(2) = 2*x(1)*x(3) + 2*x(2)*x(4);
gx(3) = 2*x(3)^2 + 2*x(4)^2 - 4/m*(x(1)^2 + x(2)^2)*x(5) - 2*g*x(2);