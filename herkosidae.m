function [ APPROX,T,H ] = herkosidae( Ab, c, s, p, x0, t0, tf, func, ...
          var, h0, delta, tol, ptol, Estat, ssc, Jopt, Nopt, eps0, beta )
%% Wrapper for core iteration of HERKosiDAE 

% Initialize approximation matrix APPROX.
n = length(x0);

% Create function handles E,f and g from given function string.
E = str2func([func,'_E']);
f = str2func([func,'_f']);
g = str2func([func,'_g']);
J = str2func([func,'_J']);

% Preallocate memory for Jacobian J. This is faster than allocating memory
% in every iteration step.
[m,~] = size(g(x0,t0,var));
Jac = zeros(m,n);

% Define variables for decomposition of leading matrix. These are
% used if option Estat equals 1.
LE = []; UE = []; PE = []; QE = []; uE = 0;

% Save dimensions for later use.
dim = [m,n];

% Use adaptive step size control if ssc is set accordingly.
if ssc == 1
    
    % Initialize approximation matrix with initial value.
    APPROX(:,1) = x0;
    T(:,1) = t0;
    H(:,1) = h0;
    t = t0;
    j = 1;
    h = h0;

    % Throw warning if tolerance for Newton iteration is larger than
    % desired accuracy.
    if tol > eps0
       warning(['Tolerance for Newton iteration is larger than ', ...
           'desired accuracy. Adjusted tol to be equal to eps0.']);
       tol = eps0;
    end
    
    % Throw warning if differentiation limit is lower than desired 
    % accuracy.
    if delta > eps0
       warning(['One should consider to lower delta. If delta is ', ...
           'larger than tol Newton iteration may take long to converge.']);
    end
    
    % Time steps
    while t < tf

        % Initialize current time step
        X0 = APPROX(:,j);
        T0 = t;

        % Half-explicit Runge-Kutta step with adaptive step control
        while 1
            
            % Take current step twice with current step size and halved
            % step size.
            [Xi, LE, UE, PE, QE, uE] = herkosidae_core(E, f, g, J, Jac, ...
                   LE, UE, PE, QE, uE, Jopt, Nopt, Ab, c, s, h, eps0, ...
                   X0, T0, Estat, delta, tol, ptol, var, dim);
            [Xi2, LE, UE, PE, QE, uE] = herkosidae_core(E, f, g, J, Jac,... 
                   LE, UE, PE, QE, uE, Jopt, Nopt, Ab, c, s, h/2, eps0, ...
                   X0, T0, Estat, delta, tol, ptol, var, dim);
            [Xi3, LE, UE, PE, QE, uE] = herkosidae_core(E, f, g, J, Jac,...
                   LE, UE, PE, QE, uE, Jopt, Nopt, Ab, c, s, h/2, eps0, ...
                   Xi2, T0 + h/2, Estat, delta, tol, ptol, var, dim);
            
            % Compare accuracy gained by halving step size.
            epsi = norm(Xi - Xi3)/(2^p - 1);
            
            % Accept current time step and increase step size.
            if epsi <= eps0
                t = t + h;  
                hnew = beta*h*(eps0/epsi)^(1/(p+1));
                if t + hnew > tf
                    h = tf - t;
                else
                    h = hnew;
                end           
                j = j + 1;
                break;
            % Deny current time step and lower step size.
            else 
                h = beta*h*(eps0/epsi)^(1/p);
            end
        end
        
        % Save current approximation and time step.
        T(:,j) = t;
        H(:,j) = h;
        APPROX(:,j) = Xi;

    end
    
% Use constant step size if ssc is not set.
else
    
    % Override eps0 with tol.
    eps0 = tol;
 
    % Initialize approximation matrix with initial value.
    Tf = round((tf - t0)/h0);
    APPROX = zeros(n,Tf+1);
    T = zeros(n,Tf+1);
    APPROX(:,1) = x0;
    T(:,1) = t0;
    H(:,1) = h0;
    
    % Time steps
    for t = 2:Tf+1

        % Initialize current time step
        X0 = APPROX(:,t-1);
        T0 = t0 + (t-2)*h0;
        
        % Half explicit Runge-Kutta step
        [Xi, LE, UE, PE, QE, uE] = herkosidae_core(E, f, g, J, Jac, ...
                LE, UE, PE, QE, uE, Jopt, Nopt, Ab, c, s, h0, eps0, ...
                X0, T0, Estat, delta, tol, ptol, var, dim);

        % Save current approximation and time step.
        T(:,t) = T0 + h0;
        H(:,t) = h0;
        APPROX(:,t) = Xi;

    end
end


end

function [ Xi, LE, UE, PE, QE, uE ] = herkosidae_core(E, f, g, J, Jac, ...
                        LE, UE, PE, QE, uE, Jopt, Nopt, Ab, c, s, h, ...
                        eps0, X0, T0, Estat, delta, tol, ptol, var, dim)
%% Core iteration of HERKosiDAE

    % Initialize Xdot. This might not be necessary. Maybe we can allocate
    % memory before.
    Xdot = zeros(dim(2),s);
     
    % Runge-Kutta steps
    for i = 1:s+1

        if i == 1
            % a) and b) for s = 1
            Xi = X0;        
        else
            % a)
            XDi = SD*(X0 + h*Xdot*Ab(i,:)');
            % b)
            if Nopt == 1
                XAi = simpledivnewton(g,Xi,XDi,J,Jac,Jopt,...
                                     q,sa,SA,T0+c(i)*h,delta,tol,var);     
            else
                XAi = divnewton(g,Xi,XDi,J,Jac,Jopt,...
                                        q,sa,SA,T0+c(i)*h,delta,tol,var);
            end
            % Combine differential and algebraic components
            Xi = XDi + XAi;
        end

        % c)
        if i < s+1
            [Xdot(:,i),FV,LE,UE,PE,QE,uE] = ...
                solveExf(E,LE,UE,PE,QE,uE,f,Xi,T0+c(i)*h,...
                                    var,dim,eps0,tol,Estat);
        end

        % Determine selectors once every time step
        if i == 1
            % Note that selectors are sparse matrices. Also, utilize
            % little trick by theoretically using SA = sa'*sa and 
            % SD = sd'*sd to remember position of components.
            [sa,SA,SD,q] = selectors(g,J,Jopt,Jac,Xi,T0+c(i)*h,...
                                            delta,ptol,FV,var,dim);
        end

    end
end