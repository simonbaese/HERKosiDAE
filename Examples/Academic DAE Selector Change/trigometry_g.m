function gx = trigometry_g( x,t,~ )
%% TRIGOMETRY_G

gx = zeros(2,1);
gx(1) = x(1)^2 + x(2)^2 - 1;
gx(2)  = 2*x(1)*x(2) - 2*sin(t)*x(2)*x(3);