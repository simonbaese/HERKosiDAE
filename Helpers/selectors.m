function [ sa,SA,SD,q ] = selectors(g,J,Jopt,Jac,x,t,delta,ptol,FV,var,dim)
%% Determine algebraic and differential selector

% Take a short cut here if the number of required algebraic components
% equals m. In this case there are no variable algebraic components
% and the rest of the components are differential components.
if sum(FV) == dim(1)
    
    % Initialize fake permutation matrix Q.
    Q = zeros(dim(2)); j = 1;
    q = zeros(1,dim(2));
    
    % Fill Q based on the information in FV.
    for i = 1:dim(2)
        if FV(i) == 1
           q(j) = i;
           Q(i,j) = 1; j = j + 1;
        else
           Q(i,dim(2)-i+2-j) = 1;
        end
    end
   
% Otherwise we will use a special LU decomposition of the Jacobian to get
% the information about the variable components.
else
    
    % Calculate Jacobian of g with func_J or by numerical approximation.
    if Jopt == 1
        Jac = J(x,t,var);
    else
        Jac = numjacobian(g,Jac,1:dim(2),x,t,delta,var);
    end
    
    % Retrieve permutation matrix Q and vector q from LU decomposition 
    % with special pivoting where LU = PAQ in matrix form.
    [Q,q] = lusp(Jac,FV,ptol);
    
end

% Split according to regular and singular part of Jacobian.
% sa is the selector for the algbraic components.
% sd is the selector for the differential components.
sa = Q(:,1:dim(1))';
sd = Q(:,dim(1)+1:end)';

% Construct inflated selectors. We use these to remember the position
% of the components. 
SA = sa'*sa;
SD = sd'*sd;

% Get column numbers for algebraic components. These will be used in the
% calculation of the Jacobian.
q = q(1:dim(1))';