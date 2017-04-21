function [ ] = run( )
%% MISC

% Add parent path.
cd('../../');
addpath(pwd);
addpath([pwd,'\helpers']);
cd('examples\Linear Circuit with one CV Loop');

%% EXAMPLE LINEAR CIRCUIT WITH ONE CV LOOP

% Set tolerances
delta = 1e-13;      % Differentiation limit
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
[Ab,c,s,p] = getRKmethod(7);

% Set ssc = 'adaptive' for adaptive step size control.
% Defaults to constant step size for any other value. In that case eps0 and
% beta are not used.
ssc = 1;
eps0 = 1e-10;        % Desired accuracy
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
func = 'circuit';

% Set var for evaluation of functions.
var = [];

% Starting value x0 = [q1,q2,e1,e2,iV].
q1 = 0; q2 = 0; 
e1 = 0; e2 = 0;
iV = -50;
x0 = [q1,q2,e1,e2,iV]';

% Initialize steps with 
% t0: initial time, 
% tf: final time and 
% h: step size.
format long;
t0 = 0;
tf = 1;
h0 = 1/1000;

% Calculate approximation.
%---------------------------------------------------------------
fprintf(['Start time: ',datestr(clock,'HH:MM:SS'),'\n']);
tic
[APPROX,T,H] = herkosidae(Ab, c, s, p, x0, t0, tf, func, var, h0, ...
                    delta, tol, ptol, Estat, ssc, Jopt, Nopt, eps0, beta);
toc
fprintf(['End time: ',datestr(clock,'HH:MM:SS'),'\n']);

% Calculate error.
%---------------------------------------------------------------
tt = length(APPROX(1,:));
ERR = zeros(1,tt);
ERRMAX = zeros(1,tt);
ERRX = zeros(1,tt);
ERRY = zeros(1,tt);
ERRV = zeros(1,tt);
ERRW = zeros(1,tt);
ERRL = zeros(1,tt);
ee = zeros(1,tt);
for t = 1:length(APPROX(1,:))
    z = T(t);
    e1 = sin(100*z);
    e2 = 100/40001*cos(100*z) + 20000/40001*sin(100*z) - 100/40001*exp(-0.5*z);
    q1 = e1 - e2;
    q2 = e2;
    ee(t) = q1;
    iV = -2000100/40001*cos(100*z) - 50001/40001*sin(100*z) + 50/40001*exp(-0.5*z);
    ERR(1,t) = norm(APPROX(:,t) - [q1,q2,e1,e2,iV]');
    ERRX(1,t) = norm(APPROX(1,t) - q1);
    ERRY(1,t) = norm(APPROX(2,t) - q2);
    ERRV(1,t) = norm(APPROX(3,t) - e1);
    ERRW(1,t) = norm(APPROX(4,t) - e2);
    ERRL(1,t) = norm(APPROX(5,t) - iV);
    ERRMAX(1,t) = max(ERR(1,1:t));
end

close all;

% Output error.
%---------------------------------------------------------------
figure('Name', 'Error linear circuit', ...
   'NumberTitle', 'off','InnerPosition',[1 1 1000 250]);
y1 = ERRMAX(1,:);
plot(T,y1,'r','LineWidth',1.3)
hold on;
y2 = ERR(1,:);
plot(T,y2,'b','LineWidth',1.2);
%title('Error linear circuit')
xlabel('time','FontWeight','bold')
ylabel('error','FontWeight','bold')
legend('max(error)','Location','northwest')
hold off;

% Plot error in components.
% ---------------------------------------------------------------
figure('Name', 'Error Linear Circuit', ...
   'NumberTitle', 'off','InnerPosition',[1 1 1000 1000]);
t = T;
x = ERRX(1,:);
subplot(5,1,1)
plot(t,x,'b','LineWidth',1.3);
ylabel('error in q_1','FontWeight','bold')
y = ERRY(1,:);
subplot(5,1,2)
plot(t,y,'b','LineWidth',1.3);
ylabel('error in q_2','FontWeight','bold')
v = ERRV(1,:);
subplot(5,1,3)
plot(t,v,'b','LineWidth',1.3);
ylabel('error in e_1','FontWeight','bold')
w = ERRW(1,:);
subplot(5,1,4)
plot(t,w,'b','LineWidth',1.3);
ylabel('error in e_2','FontWeight','bold')
z = ERRL(1,:);
subplot(5,1,5)
plot(t,z,'b','LineWidth',1.3);
ylabel('error in i_V','FontWeight','bold')
xlabel('time','FontWeight','bold')
% print('errorrk6','-dpng')

% Plot stepsize.
% ------------------------------------------------------------------
figure('Name', 'Step size Chemakzo', 'NumberTitle', 'off');
plot(T,H,'b','LineWidth',1.3);
ylabel('step size','FontWeight','bold');
xlabel('time','FontWeight','bold');

% Plot components.
% ------------------------------------------------------------------
figure('Name', 'Components Linear Circuit', 'NumberTitle', 'off');
subplot(5,1,1)
yy = APPROX(1,:);
plot(T,yy,'b','LineWidth',1.3);
title('q_1','FontWeight','bold');
subplot(5,1,2)
yy = APPROX(2,:);
plot(T,yy,'b','LineWidth',1.3);
title('q_2','FontWeight','bold');
subplot(5,1,3)
yy = APPROX(3,:);
plot(T,yy,'b','LineWidth',1.3);
title('e_1','FontWeight','bold');
subplot(5,1,4)
yy = APPROX(4,:);
plot(T,yy,'b','LineWidth',1.3);
title('e_2','FontWeight','bold');
subplot(5,1,5)
yy = APPROX(5,:);
plot(T,yy,'b','LineWidth',1.3);
title('i_V','FontWeight','bold');

% Output error and required steps.
% ---------------------------------------------------------------
format long;
fprintf('The error for the final step is %s.\n', ERR(1,end));
fprintf('Approximation required %d time steps.\n', length(T));