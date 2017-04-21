function fx = massspringchain_f( x,t,var )
%% MASSSPRINGCHAIN_F

% Retrieve C
C = var(2);

% Evaluate function
fx = zeros(7,1);

fx(1) = x(4);
fx(2) = x(5);
fx(3) = x(6);
fx(4) = -C*(x(1) - x(2)) + x(7);
fx(5) = C*(x(1) - x(2)) - C*(x(2) - x(3));
fx(6) = C*(x(2) - x(3)) + x(7);
fx(7) = x(2) - sin(t);