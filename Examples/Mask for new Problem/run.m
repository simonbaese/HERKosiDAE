function [ ] = run( )
%% MISC

% Add parent path.
cd('../../');
addpath(pwd);
addpath([pwd,'\helpers']);
cd('examples\Mask for new Problem');

%% MASK FOR NEW PROBLEM

% Set tolerances
delta = 1e-14;      % Differentiation limit
tol = 1e-7;        % Tolerance for newton iteration
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
[Ab,c,s,p] = getRKmethod(4);

% Set ssc to 1 for adaptive step size control.
% Defaults to constant step size for any other value. In that case eps0 and
% beta are not used.
ssc = 1;
eps0 = 1e-7;        % Desired accuracy
beta = 0.8;          % Safety factor for step size control

% If Jacobian is analytically known, it can be defined in func_J.m! Set
% option to 1.
Jopt = 0;

% Determine which Newton method should be used.
% Set option Nopt to 1 for simplified Newton method.
% Set option Nopt to 0 for classic Newton method.
Nopt = 1;

% If leading matrix E is time invariant, set option Estat to 1.
Estat = 1;

% Set function string
func = 'func';

% Set var for evaluation of functions
c1 = 1;
c2 = 1;
var = [c1,c2];

% Starting value.
x0 = [0,0,0,0,0]';

% Initialize steps with 
% t0: initial time, 
% tf: final time and 
% h0: initial step size.
t0 = 0;
tf = 1;
h0 = 1/100;

% Calculate approximation.
fprintf(['Start time: ',datestr(clock,'HH:MM:SS'),'\n']);
tic
[APPROX,T,H] = herkosidae(Ab, c, s, p, x0, t0, tf, func, var, h0, ...
                    delta, tol, ptol, Estat, ssc, Jopt, Nopt, eps0, beta);
toc
fprintf(['End time: ',datestr(clock,'HH:MM:SS'),'\n']);