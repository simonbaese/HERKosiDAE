function [ gx ] = circuit_g( x,t,~ )
%% CIRCUIT_G 

gx = zeros(4,1);

gx(1) = x(3) - sin(100*t);
gx(2) = x(1) - x(3) + x(4);
gx(3) = x(2) - x(4);
gx(4) = 2*x(3) + x(4) + 2*x(5) + 100*cos(100*t);