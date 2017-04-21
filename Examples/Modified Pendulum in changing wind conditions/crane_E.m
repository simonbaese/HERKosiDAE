function [ E ] = crane_E( ~,~,var )
%% CRANE_E

% Retrieve m
m = var(1);

% Evaluate function
E = zeros(5);

E(1,1) = 1;
E(2,2) = 1;
E(3,3) = m;
E(4,4) = m;