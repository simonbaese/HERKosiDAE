function Jac = numjacobian(g,Jac,q,x,t,delta,var)
%% Simple function to determine entries of Jacobian of given function.

% Evaluate given function at (x,t)
gx = g(x,t,var);

% Calculate partial derivatives
for j = q
    x(j) = x(j) + delta;   
    Jac(:,j) = (g(x,t,var) - gx) / delta;
    x(j) = x(j) - delta;   
end