function [Ab, c, s, p] = getRKmethod( run )
%% Get Runge-Kutta method

% Explicit Euler
if run == 1
    A = 0;
    c = 0;
    b = 1;
    s = 1;
    p = 1;
end

% Heun method
if run == 2
    A = [[0,0];[1,0]];
    c = [0,1];
    b = [1/2,1/2];
    s = 2;
    p = 2;
end

% Kutta's third-order method
if run == 3
    A = [[0,0,0];[1/2,0,0];[-1,2,0]];
    c = [0,1/2,1];
    b = [1/6,2/3,1/6];
    s = 3;
    p = 3;
end

% Classic fourth-order method
if run == 4
    A = [[0,0,0,0];[1/2,0,0,0];[0,1/2,0,0];[0,0,1,0]];
    c = [0,1/2,1/2,1];
    b = [1/6,1/3,1/3,1/6];
    s = 4;
    p = 4;
end

% Brasey-Hairer 3-Stage HERK
if run == 5  
    A = [[0,0,0];[1/3,0,0];[-1,2,0]];
    c = [0,1/3,1];
    b = [0,3/4,1/4];
    s = 3;
    p = 3;
end

% Brasey-Hairer 5-Stage HEM4
if run == 6
    A = [[0,0,0,0,0]; ...
        [3/10,0,0,0,0]; ...
        [(1+sqrt(6))/30,(11-4*sqrt(6))/30,0,0,0]; ...
        [(-79-31*sqrt(6))/150,(-1-4*sqrt(6))/30,(24+11*sqrt(6))/25,0,0]; ...
        [(14+5*sqrt(6))/6,(-8+7*sqrt(6))/6,(-9-7*sqrt(6))/4,(9-sqrt(6))/4,0]];
    c = [0,3/10,(4-sqrt(6))/10,(4+sqrt(6))/10,1];
    b = [0,0,(16-sqrt(6))/36,(16+sqrt(6))/36,1/9];
    s = 5; 
    p = 4;
end

% 3/8-rule fourth-order method
if run == 7
    A = [[0,0,0,0];[1/3,0,0,0];[-1/3,1,0,0];[1,-1,1,0]];
    c = [0,1/3,2/3,1];
    b = [1/8,3/8,3/8,1/8];
    s = 4;
    p = 4;   
end

% Combine A and b for Runge-Kutta iteration and
% append 1 to c for the final Runge-Kutta step
Ab = [A;b];
c = [c,1];

end

