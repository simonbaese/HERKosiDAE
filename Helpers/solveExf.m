function [ xdot,FV,L,U,P,Q,u ] = solveExf( E,LE,UE,PE,QE,uE,f,x,t,...
                                                var,dim,eps0,tol,Estat )
%% Solves E(x,t)xdot = f(x,t)

% Check if option Estat is set.
if Estat == 1 && ~isempty(LE)
    
    L = LE; U = UE; P = PE; Q = QE; u = uE;

% Otherwise we decompose E to solve the system of equations.
else
    
    % Use LU decomposition with full pivoting.
    [L,U,P,Q] = lu(sparse(E(x,t,var)));
    
    % Get rank of U.
    % Here we trust that the ranks within the model are constant. 
    % If the model is not trusted one can use u = rank(U,tol).
    u = full(sum(diag(U) ~= 0));
    
end

% Save non-zero factors of LU decomposition.
% All zero factors are related to a required algebraic component.
ft = P*f(x,t,var);

% We check if inhomogenity is fulfilled. This warning can be commented
% out when weak model formulation is tested.
if u < dim(2) && norm(ft(u+1:dim(2))) > max(100*tol,100*eps0)
    warning('Inhomogenity not fulfilled - check model');
end

% Solve main equation for xdot.
% If all entries of ft are zero fall back.
if any(ft)
    xdot = Q*[(L(1:u,1:u)*U(1:u,1:u))\ft(1:u,1);zeros(dim(2)-u,1)];
else
    xdot = zeros(dim(2),1);
end

% Get indicator of required algebraic components with permutation 
% matrix of LU decomposition.
FV = Q*[zeros(u,1);ones(dim(2)-u,1)];