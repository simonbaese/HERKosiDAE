function [ Jac ] = pendulum_J( x,~,var )
%% PENDULUM_J 

% Retrieve constants
m = var(1);
l = var(2);
g = var(3);

% Evaluate function
% Jac = [-4/m*x(1)*x(5),-4/m*x(2)*x(5)-g,2*x(3),2*x(4),-2/m*l^2];

Jac = [[2*x(1),2*x(2),0,0,0];
     [2*x(3),2*x(4),2*x(1),2*x(2),0];
     [-8/m*x(1)*x(5),-8/m*x(2)*x(5)-2*g,4*x(3),4*x(4),-4/m*l^2]];