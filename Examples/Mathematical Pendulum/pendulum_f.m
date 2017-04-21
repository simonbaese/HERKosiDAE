function fx = pendulum_f( x,~,var )
%% PENDULUM_F

% Retrieve constants
m = var(1);
l = var(2);
g = var(3);

% Evaluate function
fx = zeros(5,1);

fx(1) = x(3);
fx(2) = x(4);
fx(3) = -2*x(1)*x(5);
fx(4) = -m*g-2*x(2)*x(5);
fx(5) = x(1)^2 + x(2)^2 - l^2;