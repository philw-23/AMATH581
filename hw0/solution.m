clear;
close all;
clc

%% Exercise 1: Building a Matrix 
A = [34 45; 17 6];
A1 = A; % ANSWER

%% Exercise 2: Matrix Operations
A = [1 2; -1 1];
B = [2 0; 0 2];
C = [2 0 -3; 0 0 -1];
D = [1 2; 2 3; -1 0];
x = [1; 0];
y = [0; 1];
z = [1; 2; -1];

A2 = A + B; % ANSWER
A3 = 3*x - 4*y; % ANSWER
A4 = A * x; % ANSWER
A5 = B * (x - y); % ANSWER
A6 = D * x; % ANSWER
A7 = D * y + z; % ANSWER
A8 = A * B; % ANSWER
A9 = B * C; % ANSWER
A10 = C * D; % ANSWER

%% Exercise 3: Root Finding

% Netwon-Rhapson Method
tol = 10^-6; % Error Tolerance
x_0 = -3; % Initial Condition
f = @(x) -x - cos(x); % Function
f_prime = @(x) -1 + sin(x); % Derivative
f_vec_n = []; % Solution vector (Don't include initial condition)
% NOTE: We start at 0 because we are storing x_values
% The first value to store is x_0, thus there it no iteration to generate
% it. Therefore the second value (x_next) is essentially the first
% iteration
iter_n = 0; % Starting iteration
max_iter = 1000; % Max iterations

x_next = x_0 - f(x_0) / f_prime(x_0); % First iteration
err = abs(x_next - x_0); % Initial error
% NOTE: Want to store x, not x_next
% x_next is for error
x = x_0; % Convenience Variable
f_vec_n(end + 1) = x; % Append

while (err > tol) && (iter_n <= max_iter)
    iter_n = iter_n + 1;
    x = x_next; % Update
    x_next = x - f(x) / f_prime(x); % Forward iteration
    f_vec_n(end + 1) = x;
    err = abs(x_next - x);
end

% Bisection Method
x_a = -3; % First endpoint
x_b = 1; % Second endpoint
iter_b = 1; % Bisection method iterations
c = (x_a + x_b) / 2; % First midpoint
f_mid = f(c); % Function  evaluated at midpoint
f_vec_b = [c]; % Solution Vector

while (f_mid ~= 0) && ((x_b - x_a)/2 > tol) && (iter_b <= max_iter)
    if sign(f(c)) == sign(f(x_a))
        x_a = c;
    else
        x_b = c;
    end
    c = (x_a + x_b)/2;
    f_mid = f(c);
    f_vec_b(end + 1) = c;
    iter_b = iter_b + 1; % increment
end

A11 = f_vec_n'; % ANSWER
A12 = f_vec_b'; % ANSWER
A13 = [iter_n iter_b]; % ANSWER
