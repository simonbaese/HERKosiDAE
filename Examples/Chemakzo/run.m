function [ ] = run( )
%% MISC

% Add parent path.
cd('../../');
addpath(pwd);
addpath([pwd,'\helpers']);
cd('examples\Chemakzo');

%% EXAMPLE CHEMAKZO

% Set tolerances
delta = 1e-14;      % Differentiation limit
tol = 1e-6;        % Tolerance for newton iteration
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
[Ab,c,s,p] = getRKmethod(6);

% Set ssc = 'adaptive' for adaptive step size control.
% Defaults to constant step size for any other value. In that case eps0 and
% beta are not used.
ssc = 1;
eps0 = 1e-6;        % Desired accuracy
beta = 0.78;          % Safety factor for step size control

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
func = 'chemakzo';

% Set var for evaluation of functions.
k1 = 18.7;
k2 = 0.58;
k3 = 0.09;
k4 = 0.42;
K = 34.4;
klA = 3.3;
Ks = 115.83;
pCO2 = 0.9;
H = 737;
var = [k1,k2,k3,k4,K,klA,Ks,pCO2,H];

% Starting value x0 = [y1,y2,y3,y4,y5,y6]
y1 = 0.444; y2 = 0.00123; y3 = 0;
y4 = 0.007; y5 = 0; y6 = Ks*y1*y4;
x0 = [y1,y2,y3,y4,y5,y6]';

% Initialize steps with 
% t0: initial time, 
% tf: final time and 
% h: step size.
t0 = 0;
tf = 180;
h0 = 1/100;

% Calculate approximation.
fprintf(['Start time: ',datestr(clock,'HH:MM:SS'),'\n']);
tic
[APPROX,T,H] = herkosidae(Ab, c, s, p, x0, t0, tf, func, var, h0, ...
                    delta, tol, ptol, Estat, ssc, Jopt, Nopt, eps0, beta);
toc
fprintf(['End time: ',datestr(clock,'HH:MM:SS'),'\n']);


% Output error to reference solution and required steps.
% ------------------------------------------------------------------
format long;
SOL = [0.1150794920661702,...
       0.0012038314715677, ...
       0.1611562887407974, ...
       0.0003656156421249, ...
       0.0170801088526440, ...
       0.0048735313103074]';

fprintf('The error for the final step is %s.\n', ...
    norm(APPROX(:,end)-SOL));
fprintf('Approximation required %d time steps.\n', length(T));

% Plot stepsize.
% ------------------------------------------------------------------
figure('Name', 'Step size Chemakzo', 'NumberTitle', 'off');
plot(T,H,'b','LineWidth',1.3);
ylabel('step size','FontWeight','bold');
xlabel('time','FontWeight','bold');

% Plot components.
% ------------------------------------------------------------------
figure('Name', 'Components Chemakzo', 'NumberTitle', 'off');
subplot(3,2,1)
yy = APPROX(1,:);
plot(T,yy,'b','LineWidth',1.3);
title('y_1','FontWeight','bold');
subplot(3,2,2)
yy = APPROX(2,:);
plot(T,yy,'b','LineWidth',1.3);
title('y_2','FontWeight','bold');
subplot(3,2,3)
yy = APPROX(3,:);
plot(T,yy,'b','LineWidth',1.3);
title('y_3','FontWeight','bold');
subplot(3,2,4)
yy = APPROX(4,:);
plot(T,yy,'b','LineWidth',1.3);
title('y_4','FontWeight','bold');
subplot(3,2,5)
yy = APPROX(5,:);
plot(T,yy,'b','LineWidth',1.3);
title('y_5','FontWeight','bold');
subplot(3,2,6)
yy = APPROX(6,:);
plot(T,yy,'b','LineWidth',1.3);
title('y_6','FontWeight','bold');