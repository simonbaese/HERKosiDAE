function [ fx ] = circuit_f( x,t,~ )
%% CIRCUIT_F

fx = zeros(5,1);
fx(1) = x(3) + x(5);
fx(2) = x(4);
fx(3) = x(3) - sin(100*t);
fx(4) = x(1) - x(3) + x(4);
fx(5) = x(2) - x(4);