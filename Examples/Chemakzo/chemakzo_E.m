function [ E ] = chemakzo_E( ~,~,~ )
%% CHEMAKZO_E

E = eye(6);
E(6,6) = 0;