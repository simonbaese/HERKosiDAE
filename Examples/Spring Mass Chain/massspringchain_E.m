function [ E ] = massspringchain_E( ~,~,var )
%% MASSSPRINGCHAIN_E

% Retrieve m
m = var(1);

% Evaluate function
E = zeros(7);

E(1,1) = 1;
E(2,2) = 1;
E(3,3) = 1;
E(4,4) = m;
E(5,5) = m;
E(6,6) = m;