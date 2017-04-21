function [ Ex ] = func_E( x,t,var )
%% FUNC_E

% Retrieve constants
c1 = var(1);
% c2 = var(2);

Ex = zeros(5);

% Evaluate function
Ex(1,1) = 1;
Ex(2,2) = 1;
Ex(3,3) = c1;
Ex(4,4) = c1;