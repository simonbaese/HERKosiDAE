function XAk = simpledivnewton( g,X,XD,J,Jac,Jopt,...
                                q,sa,SA,t,delta,tol,var )
%% Simple Devided Newton method

% Initialise selected algebraic part
XAk = SA*X;

% Test if any components are algebraic.
if ~isempty(q)

    % Get selected Jacobian with func_J or numerical approximation. 
    % Here we add the differential part from step a).
    if Jopt == 1
        JAk = J(XD+XAk,t,var)*sa';
    else
        JAk = numjacobian(g,Jac,q,XD+XAk,t,delta,var)*sa';
    end
    
    % Decompose Jacobian to quickly calculate inverse.
    [Q,R,P] = qr(JAk);

    % Test condition of matrix.
    COND = R(1,1)/R(end,end);
    if abs(COND) > 1e5
        warning('Jacobian ill conditioned!');
    end

    % Get factor, inverse only once.
    C = P/R;

    % Simplified Newton approximation.
    for k = 1:1000

        % Keep differential part and evaluate function.
        Gk = g(XD+XAk,t,var);

        % Check convergence.
        if abs(sum(Gk)) < tol
            break;
        elseif k == 100
            % Try to be graceful
            if abs(sum(Gk)) < 10*tol
                warning('Newton method reached only precision of 10*tol!');
                break;
            else
                error('Newton iteration does not converge!');
            end
        end

        % Progress algebraic part.
        % This value is returned once loop breaks.
        XAk = XAk - sa'*C*Q'*Gk;

    end
end