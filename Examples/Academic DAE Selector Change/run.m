function [] = run()
%% MISC

% Add parent path.
cd('../../');
addpath(pwd);
addpath([pwd,'\helpers']);
cd('examples\Academic DAE Selector Change');

%% EXAMPLE TRIGONOMETRY

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
[Ab,c,s,p] = getRKmethod(3);

% Set ssc = 'adaptive' for adaptive step size control.
% Defaults to constant step size for any other value. In that case eps0 and
% beta are not used.
ssc = 1;
eps0 = 1e-12;        % Desired accuracy
beta = 0.9;          % Safety factor for step size control

% When Jacobian is analytically known, it can be defined in func_J.m! 
% Set option jopt to 1 in that case.
% Set option jopt to 0 otherwise.
Jopt = 0;

% Determine which Newton method should be used.
% Set option Nopt to 1 for simplified Newton method.
% Set option Nopt to 0 for classic Newton method.
Nopt = 0;

% If leading matrix E is time invariant, set option Estat to 1.
Estat = 1;

% Set function string.
func = 'trigometry';

% Set var for evaluation of functions.
var = [];

% Set initial value x0.
x0 = [sin(pi/8),cos(pi/8),1]';

% Initialize steps with 
% t0: initial time, 
% tf: final time and 
% h: step size.
% Note that Jacobian in step b) will be singular for k*pi/2 with k integer.
t0 = pi/8;
tf = 3*pi/8;
h0 = (3*pi/8-pi/8)/100;

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
SOL = [sin(tf),cos(tf),1]';
fprintf('The error for the final step is %s.\n', ...
    norm(APPROX(:,end)-SOL));
fprintf('Approximation required %d time steps.\n', length(T));

% Output plot.
%---------------------------------------------------------------
figure('Name', 'Academic Example', 'NumberTitle', 'off');
y1 = APPROX(1,:);
y2 = APPROX(2,:);
y3 = APPROX(3,:);
plot(T,y1,T,y2,T,y3,'LineWidth',1.3);