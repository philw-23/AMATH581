clear;
close;
clc;

%% Exercise 1a - Euler's Method

dydt = @(t, y) -3 * y * sin(t); % dy/dt
y_exact = @(t) (pi * exp(3 * (cos(t) - 1))) / sqrt(2); % Exact sol
y_0 = pi / sqrt(2); % initial condition
y_n = y_0; % Convenience variable for starting
j = -2:-1:-8;
delta_t = 2.^j; % Delta t values to consider
t_bounds = [0 5];
l_dt = length(delta_t); % Helper variable
err_vec = zeros(l_dt, 2); % initialize error vector
err_vec(:, 1) = delta_t;

for i = 1:l_dt
    dt = delta_t(i); % Delta t
    t = [0:dt:5]; % t vector
    y_true = y_exact(t); % True results
    y_pred = []; % Prediction
    for j = 1:length(t)
        y_new = y_n + dt * dydt(t(j), y_n); % Iterate
        y_pred(end + 1) = y_n; % add to results
        y_n = y_new; % Update
    end
    y_n = y_0; % reset
    e_dt = mean(abs(y_true - y_pred)); % error for this iteration
    err_vec(i, 2) = e_dt;
end

log_err = log(err_vec); % Take the log
p = polyfit(log_err(:, 1), log_err(:, 2), 1); % linear_fit

f = figure('visible', 'off'); % figure object
plot(log_err(:, 1), log_err(:, 2)) % Plot for reference
xlabel('Log of delta t')
ylabel('Log of Error')
saveas(f, 'Euler_Method_Error.png')

A1 = y_pred; % ANSWER - Last dt values
A2 = err_vec(:, 2)'; % ANSWER - Error values
A3 = p(1); % ANSWER - Slope of polyfit line

%% Exercise 1b - Heun's method

for i = 1:l_dt
    dt = delta_t(i);
    t = [0:dt:5];
    y_true = y_exact(t);
    y_pred = []; % prediction
    for j = 1:length(t)
        y_new = y_n + (dt / 2) * (dydt(t(j), y_n) + ...
            dydt(t(j) + dt, y_n + dt * dydt(t(j), y_n))); % iterate
        y_pred(end + 1) = y_n;
        y_n = y_new;
    end
    y_n = y_0; % reset
    e_dt = mean(abs(y_true - y_pred));
    err_vec(i, 2) = e_dt;
end

log_err = log(err_vec); % Take the log
p = polyfit(log_err(:, 1), log_err(:, 2), 1); % linear_fit

f = figure('visible', 'off'); % figure object
plot(log_err(:, 1), log_err(:, 2)) % Plot for reference
xlabel('Log of delta t')
ylabel('Log of Error')
saveas(f, 'Heun_Method_Error.png')

A4 = y_pred; % ANSWER - Last dt values
A5 = err_vec(:, 2)'; % ANSWER - Error values
A6 = p(1); % ANSWER - Slope of polyfit line

%% Exercise 2a - van der Pol Oscillator

% Define parameters
eps = [0.1 1 20]; % eps parameters
y_0 = sqrt(3); % function initial condition
dy_0 = 1; % derivative initial condition
y_sols = [];

% Write problem as series of first order ODEs
% let y_1 = y'
% y_1' = y'' =  -eps * (y^2 - 1) * y' - y
%            = -eps * (y^2 - 1) * y_1 - y
% Define system in odefun format
% y(2) represents y_1 (y'), y(1) is y
% Thus will pass initial conditions as [y_0; dy_0]
for e = eps % Solve for all epsilon values
    vdp_ode = @(t, y) [y(2); -e * (y(1)^2 - 1) * y(2) - y(1)]; 
    [t_sol, y_sol] = ode45(vdp_ode, [0:0.5:32], [y_0; dy_0]);
    y_sols(:, end + 1) = y_sol(:, 1);
end

A7 = y_sols; % ANSWER - results for y(t) for each epsilon

%% Exercise 2b - van de Pol Oscillator (variable step size)

t_span = [0 32]; % No step size specified
y_0 = 2; % initial condition
dy_0 = pi^2; % initial condition
eps = 1; % epsilon
tol_powers = -4:-1:-10; % power values for tolerance vectors
tols = 10.^tol_powers; % tolerance vectors
% define odefunc once as we are not changing parameters
vdp_ode = @(t, y) [y(2); -eps * (y(1)^2 - 1) * y(2) - y(1)];
ode_funcs = {@ode45 @ode23 @ode113}; % funcs to try
plot_titles = {'ode45' 'ode23' 'ode113'}; % Plot names
p_sols = []; % Polyfit solutions
dt_sols = []; % dt solution vector

for i = 1:length(ode_funcs)
    [p, dt, f] = runVDPOde(ode_funcs{i}, vdp_ode, t_span, ...
                            tols, [y_0 dy_0]);
    p_sols(end + 1) = p(1); % Slope of line
    dt_sols(:, end + 1) = dt; % dt_vec
    saveas(f, strcat('vdp_logplot_', plot_titles{i}, '.png'))
end

A8 = p_sols(1); % ANSWER - polyfit of ode45 tol/dt
A9 = p_sols(2); % ANSWER - polyfit of ode23 tol/dt
A10 = p_sols(3); % ANSWER - polyfit of ode113 tol/dt


