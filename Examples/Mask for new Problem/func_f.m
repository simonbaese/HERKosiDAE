function fx = func_f( x,t,var )
%% FUNC_F

% Retrieve variables
c1 = var(1);
% c2 = var(2);

% Evaluate function
fx = zeros(5,1);

fx(1) = 1;
fx(2) = 1;
fx(3) = c1;
fx(4) = c1;