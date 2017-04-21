function fx = trigometry_f( x,t,~ )
%% TRIGOMETRY_F

fx = zeros(3,1);
fx(1) = x(2);
fx(2) = -sin(t)*x(3);
fx(3) = x(1)^2 + x(2)^2 - 1;