function fx = crane_f( x,t,var )
%% CRANE_F

% Retrieve constants
m = var(1);
g = var(2);

alpha = 1 - 0.04*t;
l = 0.5 + 0.03*(50-t);

% Evaluate function
fx = zeros(3,1);

fx(1) = x(3);
fx(2) = x(4);
fx(3) = -2*x(1)*x(5) + alpha*m;
fx(4) = -m*g-2*x(2)*x(5);
fx(5) = x(1)^2 + x(2)^2 - l^2;