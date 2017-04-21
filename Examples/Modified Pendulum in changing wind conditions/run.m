function [] = run()
%% MISC

% Add parent path.
cd('../../');
addpath(pwd);
addpath([pwd,'\helpers']);
cd('examples\Modified Pendulum in changing wind conditions');

%% EXAMPLE MODIFIED PENDULUM IN CHANGING WIND CONDITIONS

% Set tolerances.
delta = 1e-12;      % Differentiation limit
tol = 1e-10;        % Tolerance for newton iteration
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
ssc = 0;
eps0 = 1e-6;        % Desired accuracy
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

% Set var for evaluation of functions.
% Length l and wind alpha will be calculated within the functions.
m = 2;          %mass
l = 2;          %this is l at t = 0
g = 9.78;       %gravitational force
var = [m,g];

% Set function string.
func = 'crane';

% Starting value x0 = [x,y,v,w,lambda] with
% free initial values
x = 0;
v = 0;

% and dependent initial values.
if x ~= abs(l)
    y = -sqrt(l^2 - x^2);
    w = -x*v/y;
else
    y = 0;
    w = 0; %free
end
lambda = (v^2 + w^2 - g*y)*m/2/l^2;
x0 = [x,y,v,w,lambda]';

% Initialize steps with 
% t0: initial time, 
% tf: final time and 
% h: step size.
% Note that if these parameters change, l(t) and alpha(t) should be 
% adjusted accordingly in crane_g and crane_f.
t0 = 0;
tf = 50;
h0 = 1/50;

% Calculate approximation.
%---------------------------------------------------------------
fprintf(['Start time: ',datestr(clock,'HH:MM:SS'),'\n']);
tic
[APPROX,T,~] = herkosidae(Ab, c, s, p, x0, t0, tf, func, var, h0, ...
                    delta, tol, ptol, Estat, ssc, Jopt, Nopt, eps0, beta);
toc
fprintf(['End time: ',datestr(clock,'HH:MM:SS'),'\n']);

% Output time steps.
%---------------------------------------------------------------
fprintf('Approximation required %d time steps.\n', length(T));

% Output plot.
%---------------------------------------------------------------
figure('Name', 'Crane Pendulum', ...
   'NumberTitle', 'off','InnerPosition',[1 1 500 1000]);
x = APPROX(1,:);
y = APPROX(2,:);
plot(x,y,'b','LineWidth',1.3,'MarkerEdgeColor','b')
axis([-0.42 0.42 -2 -0.4]);
% title('Crane Pendulum position')
xlabel('x','FontWeight','bold')
ylabel('y','FontWeight','bold')