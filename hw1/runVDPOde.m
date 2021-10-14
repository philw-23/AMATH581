% run VDP ode function
% returns: 
%           p = error of log(delta t) vs log(error)
%           dt = vector of mean step sizes
%           f = figure object plotting log(delta t) vs log(error)
% inputs: 
%           ode_fun = sepcific ode solver to use (ode45, ode113, etc.)
%           ode = indivual ode function to solve
%           t_span = time span to iterate through
%           tols = tolerance vectors
%           init_cond = initial conditions for ode
function [p, dt, f] = runVDPOde(ode_fun, ode, t_span, tols, init_cond)
    
    dt = []; % vector for dt vals
    y_0 = init_cond(1); % y_0 initial condition
    dy_0 = init_cond(2); % dy_0 initial condition

    for t = tols % iterate though 
        options = odeset('AbsTol', t, 'RelTol', t); % tolerance options
        % run ode45 with tolerance options
        [t_sol, y_sol] = ode_fun(ode, t_span, [y_0; dy_0], options);
        mean_step = mean(diff(t_sol)); % mean delta t step
        dt(end + 1) = mean_step; % mean delta t
    end
    
    % Plot
    log_vals = log([dt' tols']); % take log
    f = figure('visible', 'off');
    plot(log_vals(:, 1), log_vals(:, 2));
    xlabel('Log of delta t');
    ylabel('Log of TOL')
   
    p = polyfit(log_vals(:, 1), log_vals(:, 2), 1); % linear fit

end