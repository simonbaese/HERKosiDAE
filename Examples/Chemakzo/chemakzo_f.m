function [ fy ] = chemakzo_f( y,~,var )
%% CHEMAKZO_F

% Retrieve variables
k1 = var(1);
k2 = var(2);
k3 = var(3);
k4 = var(4);
K = var(5);
klA = var(6);
Ks = var(7);
pCO2 = var(8);
H = var(9);

% Catch error
if y(2) < 0
    error('Component y(2) dropped below zero. Can not calculate f(y)!');
end

% Predefine
r1 = k1 * y(1)^4 * sqrt(y(2));
r2 = k2 * y(3) * y(4);
r3 = k2/K * y(1) * y(5);
r4 = k3 * y(1) * y(4)^2;
r5 = k4 * y(6)^2 * sqrt(y(2));
Fin = klA * (pCO2/H - y(2));

% Evaluate function
fy = zeros(6,1);

fy(1) = -2*r1 + r2 - r3 - r4;
fy(2) = -0.5*r1 - r4 - 0.5*r5 + Fin;
fy(3) = r1 - r2 + r3;
fy(4) = -r2 + r3 - 2*r4;
fy(5) = r2 - r3 + r5;
fy(6) = Ks*y(1)*y(4) - y(6);