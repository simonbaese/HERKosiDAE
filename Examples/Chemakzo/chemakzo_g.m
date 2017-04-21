function [ gy ] = chemakzo_g( y,~,var )
%% CHEMAKZO_G

% Retrieve variable
Ks = var(7);

% Evaluate function
gy = Ks*y(1)*y(4) - y(6);