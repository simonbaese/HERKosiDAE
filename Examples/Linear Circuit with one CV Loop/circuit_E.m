function [ E ] = circuit_E( ~,~,~ )
%% CIRCUIT_E

E = zeros(5);
E(1,1) = -1;
E(2,1) = 1;
E(2,2) = -1;