function [] = run()
%% MISC

% Add parent path.
cd('../../');
addpath(pwd);
addpath([pwd,'\helpers']);
cd('examples\Academic Pure Algebraic Equation');

%% EXAMPLE SIMPLA

% Set tolerances.
delta = 1e-15;      % Differentiation limit
tol = 1e-15;        % Tolerance for newton iteration
ptol = 1e-15;       % Tolerance for pivots in lusp

% Get Runge-Kutta method with parameters
% 1: Forward Euler method,
% 2: Heun's method,
% 3: Kutta's third-order method,
% 4: Classic fourth-order method,
% 5: Brasey-Hairer 3-Stage HERK method,
% 6: Brasey-Hairer 5-Stage HEM4 method or
% 7: 3/8-rule fourth-order method.
% Returns Butcher-Tableau Ab, c, stages s and convergence order p.
[Ab,c,s,p] = getRKmethod(2);

% Set ssc = 'adaptive' for adaptive step size control.
% Defaults to constant step size for any other value. In that case eps0 and
% beta are not used.
ssc = 0;
eps0 = 1e-15;        % Desired accuracy
beta = 0.9;          % Safety factor for step size control

% If Jacobian is analytically known, it can be defined in func_J.m! Set
% option to 1.
Jopt = 0;

% Determine which Newton method should be used.
% Set option Nopt to 1 for simplified Newton method.
% Set option Nopt to 0 for classic Newton method.
Nopt = 1;

% If leading matrix E is time invariant, set option Estat to 1.
Estat = 1;

% Set function string.
func = 'simpla';

% Set var for evaluation of functions.
var = [];

% Set initial value x0.
x0 = 1;

% Initialize steps with 
% t0: initial time, 
% tf: final time and 
% h: step size.
t0 = 0;
tf = 1;
h0 = 1/100;

% Calculate approximation.
%---------------------------------------------------------------
fprintf(['Start time: ',datestr(clock,'HH:MM:SS'),'\n']);
tic
[APPROX,T,~] = herkosidae(Ab, c, s, p, x0, t0, tf, func, var, h0, ...
                    delta, tol, ptol, Estat, ssc, Jopt, Nopt, eps0, beta);
toc
fprintf(['End time: ',datestr(clock,'HH:MM:SS'),'\n']);

% Output error and required steps.
%---------------------------------------------------------------
format long;
SOL = exp(tf);
fprintf('The error for the final step is %s.\n', ...
    abs(APPROX(:,end)-SOL));
fprintf('Approximation required %d time steps.\n', length(T));

% Output plot.
%---------------------------------------------------------------
figure('Name', 'Academic Example', 'NumberTitle', 'off');
y = APPROX(1,:);
plot(T,y,'b','LineWidth',1.3);