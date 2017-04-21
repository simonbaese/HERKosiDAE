function [] = run()
%% MISC

% Add parent path.
cd('../../');
addpath(pwd);
addpath([pwd,'\helpers']);
cd('examples\Spring Mass Chain');

%% EXAMPLE MASS SPRING CHAIN

% Set tolerances.
delta = 1e-14;      % Differentiation limit
tol = 1e-4;        % Tolerance for newton iteration
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

% Set ssc = 'adaptive' for adaptive step size control.
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

% Set function string.
func = 'massspringchain';

% Set var for evaluation of functions.
m = 1;          %mass
C = 1/6;        %spring stiffness
var = [m,C];

% Starting value x0 = [p1,p2,p3,v1,v2,v3,F].
p1 = 0; p2 = 0; p3 = 0;
v1 = -2; v2 = 1; v3 = -2;
F = 0;
x0 = [p1,p2,p3,v1,v2,v3,F]';

% Initialize steps with 
% t0: initial time, 
% tf: final time and 
% h: step size.
t0 = 0;
tf = 400;
h0 = 1/1000;

% Calculate approximation.
%---------------------------------------------------------------
fprintf(['Start time: ',datestr(clock,'HH:MM:SS'),'\n']);
tic
[APPROX,T,~] = herkosidae(Ab, c, s, p, x0, t0, tf, func, var, h0, ...
                    delta, tol, ptol, Estat, ssc, Jopt, Nopt, eps0, beta);
toc
fprintf(['End time: ',datestr(clock,'HH:MM:SS'),'\n']);

% Calculate error.
%---------------------------------------------------------------
ERR = zeros(1,length(APPROX(1,:)));
ERRMAX = zeros(1,length(APPROX(1,:)));
ERR(1,1) = norm(x0 - [-2*sin(t0),sin(t0),-2*sin(t0),...
                      -2*cos(t0),cos(t0),-2*cos(t0),3/2*sin(t0)]');
steps = T;
for t = 1:length(APPROX(1,:))
    z = steps(t);
    ERR(1,t) = norm(APPROX(:,t) - ...
                [-2*sin(z),sin(z),-2*sin(z), ...
                 -2*cos(z),cos(z),-2*cos(z),3/2*sin(z)]');
    ERRMAX(1,t) = max(ERR(1,1:t));
end

% Output error plot.
%---------------------------------------------------------------
figure('Name', 'Error mass-spring chain', ...
   'NumberTitle', 'off','InnerPosition',[0 0 1000 250]);
x = T;
y1 = ERRMAX(1,:);
plot(x,y1,'r','LineWidth',1.3)
hold on;
y2 = ERR(1,:);
plot(x,y2,'b','LineWidth',1.2);
% title('Error mass-spring chain')
ylabel('error','FontWeight','bold')
legend('max(error)','Location','northwest')
hold off;

% Output plot.
%---------------------------------------------------------------
figure('Name', 'Mass-Spring Chain', ...
   'NumberTitle', 'off','InnerPosition',[0 0 1000 1000]);

x = T;
y1 = APPROX(1,:);
subplot(7,1,1)
plot(x,y1,'b','LineWidth',1.3,'MarkerEdgeColor','b')
axis([t0 tf -2.2 2.2]);
ylabel('p_1','FontWeight','bold')
% title('Mass-spring chain')

y2 = APPROX(2,:);
subplot(7,1,2)
plot(x,y2,'b','LineWidth',1.3,'MarkerEdgeColor','b')
axis([t0 tf -1.1 1.1]);
ylabel('p_2','FontWeight','bold')

y3 = APPROX(3,:);
subplot(7,1,3)
plot(x,y3,'b','LineWidth',1.3,'MarkerEdgeColor','b')
axis([t0 tf -2.2 2.2]);
ylabel('p_3','FontWeight','bold')

y4 = APPROX(4,:);
subplot(7,1,4)
plot(x,y4,'b','LineWidth',1.3,'MarkerEdgeColor','b')
axis([t0 tf -2.2 2.2]);
ylabel('v_1','FontWeight','bold')

y5 = APPROX(5,:);
subplot(7,1,5)
plot(x,y5,'b','LineWidth',1.3,'MarkerEdgeColor','b')
axis([t0 tf -1.1 1.1]);
ylabel('v_2','FontWeight','bold')

y6 = APPROX(6,:);
subplot(7,1,6)
plot(x,y6,'b','LineWidth',1.3,'MarkerEdgeColor','b')
axis([t0 tf -2.2 2.2]);
ylabel('v_3','FontWeight','bold')

y7 = APPROX(7,:);
subplot(7,1,7)
plot(x,y7,'b','LineWidth',1.3,'MarkerEdgeColor','b')
axis([t0 tf -1.6 1.6]);
xlabel('time','FontWeight','bold')
ylabel('F','FontWeight','bold')

% Output error and required steps.
%---------------------------------------------------------------
format long;
fprintf('The maximal error for the given time frame is %s.\n', max(ERR(1,:)));
fprintf('The error for the final step is %s.\n', ERR(1,end));
fprintf('Approximation required %d time steps.\n', length(T));