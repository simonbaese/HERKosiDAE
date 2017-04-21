function [ Ex ] = pendulum_E( ~,~,var )
%%  PENDULUM_E

% Retrieve m
m = var(1);

% Evaluate function
Ex = zeros(5);

Ex(1,1) = 1;
Ex(2,2) = 1;
Ex(3,3) = m;
Ex(4,4) = m;