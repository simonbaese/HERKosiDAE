function XAk = divnewton( g,Xi,XD,J,Jac,Jopt,q,sa,SA,t,delta,tol,var )
%% Devided Newton method

%Initialise selected parts
XAk = SA*Xi;

% Test if any components are algebraic.
if ~isempty(q)
    
    % Newton iteration
    for k = 1:1000

        % Get selected Jacobian with func_J or numerical approximation.
        % Here we add the differential part from step a).
        if Jopt == 1
            JAk = J(XD+XAk,t,var)*sa';
        else
            JAk = numjacobian(g,Jac,q,XD+XAk,t,delta,var)*sa';
        end

        % Evaluate function.
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
        XAk = XAk - sa'*(JAk\Gk);
    end
    
end