clear;
close all;
clc;
plotting = true; % control if we want to show plots

%% Exercise 1 - Harmonic Trapping Potential

% initialize parameters
% x1 = phi, x2 = phi'
% thus x1' = phi' = x2, x2 = phi''
% phi'' = [Kx^2 - eps] * phi = [Kx^2 - eps]*x1
% By problem definition we know eps > 0
maxiter = 1000; % How many iterations to consider
modes = 5; % how many solutions we are looking for
eps_0 = 0.25; % starting condition for epsilon
tol = 10^(-4); % solution tolerance
L = 4; % what we consider to be boundary
A = 1; % derivative value for shooting method at x=-L
K = 1; % k parameter
dx = 0.1; % step size in x
x_0 = [A A*sqrt(K*(L^2) - eps)]; % initial conditions: ADD NOTES ON DERIVING
x_span = [-L:dx:L]; % xspan for vectors
err = 1; % initialize error
eps_sols = zeros(5, 1); % Solutions for epsilon values
eig_sols = zeros(length(x_span), modes); % for eigenvectors
plot_keys = {}; % for labels

% Iterate for solutions
if plotting
    f = figure("Name", "Shooting Method");
    hold on
end

eps_start = eps_0; % initialize for mode
for n=1:modes % modes loop
    eps = eps_start; % intialize epsilon
    eps_step = eps_0 / 50; % default step size
    x_start = x_0; % reset initial condition
    for i=1:maxiter
        % Note: x_span is our tspan for ode45
        % x_0 is our initial condition
        [x, y] = ode45(@(x, y) oscillatorShoot(x, y, K, eps), ...
                                                x_span, x_start); %, options);
        % Check for convergence:
        % Condition: phi'(L) + sqrt(K*L^2 - eps) = 0 
        % Add negative one power to account for switches
        % This makes sure we don't go back to our first mode after
        % success
        check = y(end, 2) + sqrt(K*(L^2) - eps)*y(end, 1);
        if abs(check) < tol
            break % End loop
        elseif (-1)^(n + 1) * check > 0 
            eps = eps + eps_step; % we are too low, need to be higher
        else % we are too high, need to decrease
            eps = eps - eps_step; % Update epsilon
            eps_step = eps_step / 2; % decrease eps_step
        end
        x_start = [A A*sqrt(K*(L^2) - eps)]; % update initial condition
    end
    eps_sols(n, 1) = eps; % save solution
    % Was using eps_start instead of eps for below, which was why 
    % solutions were diverging. We want to increment off of our previous epsilon
    % solution for the next starting point
    eps_start = eps + 0.1; % pick new starting value for epsilon
    y_norm = abs(y(:, 1) / sqrt(trapz(x, y(:, 1).*(y(:, 1))))); % normalize
    eig_sols(:, n) = y_norm; % append solution
    if plotting
        plot_keys{end + 1} = strcat("mode=", string(n)); % keys for plot
        plot(x, y_norm);
    end
end

if plotting
    hold off
    legend(plot_keys)
end

% ANSWERS
A1 = eig_sols(:, 1); % first eigenvector
A2 = eig_sols(:, 2); % second eigenvector
A3 = eig_sols(:, 3); % third eigenvector
A4 = eig_sols(:, 4); % fourth eigenvector
A5 = eig_sols(:, 5); % fifth eigenvector
A6 = eps_sols; % eigenvalues

%% Exercise 2 - Direct Method

% Resusing parameters from 1
% Hint 1 - define matrix without first and last points
%{ 
   We will use a central difference theorem for these points
    phi'' - [Kx^2 - eps]*phi = 0
    phi'' - Kx^2 * phi = -eps * phi
   
   Can esimate second derivative
    phi'' ~= (phi(x_i+1) - 2phi(x_i) + phi(x_i-1))/dx^2 - ...
                                    Kx_i^2 * phi(x_i) = -eps*phi(x_i)

    -> -phi(x_i-1) + (2 + dx^2*K*x_i^2)*phi(x_i) - phi(x_i+1) = eps*phi(x_i)
   
   Using this equation, can build a diagonal matrix using only interior
   points (no x_0, x_n): -1 on diagonal -1 and diagonal +1, and then the
   (2 + ...) on the main diagonal. 

   Note this doesn't include the boundary points, we will need to overwrite 
   these points later with other values. This is because our vector for phi
   is only considering values x_2 through x_n-1, so we no don't have access
   to the point x_0 or x_n for the center difference scheme. We will need
   to derive substitutes for these values
%}
x_sub = x_span(2:end-1);
A_sub = diag(-1*ones(length(x_sub)-1, 1), -1) + ...
        diag(2 + dx^2 * K * x_sub.^2) + ...
        diag(-1*ones(length(x_sub)-1, 1), 1);

% Now we need to handle boundary condition
%{
  To substitute for x_0 at our boundary, we start with a forward difference
  scheme for x_0
    phi'(x_0) = sqrt(K*x_0^2 - eps)*phi(x_0) = ...
    (-3*phi(x_0) + 4*phi(x_1) - phi(x_2)) / 2*dx
  
  Per the problem definition, we should assume dx is small such that
  when multiplying through the term will go to 0. Thus:
    3*phi(x_0) = 4*phi(x_1) - phi(x_2)

  Using this, we can substitute into our central difference scheme for phi(x_1) 
  a value for phi(x_0) because the first row of our matrix corresponds to x_1 
  and we do not have an x_0 term. This allows us to get this value in 
  terms of phi(x_1) and phi(x_2) to plug in to our A matrix
    -4/3*phi(x_1) + 1/3*phi(x_2)^2 + 2phi(x_1) + K*dx^2*x_1^2*phi(x_1)
  -> (2/3 + K*dx^2*x_1^2)*phi(x_1) - 2/3*phi(x_2) = eps*phi(x_1)

  Similarly, for the final row using a backwards difference scheme we will
  get the following:
    -2/3*phi(x_n-2) + (2/3 + K*dx^2*(x_n-1)^2)*phi(x_n-1) = eps*phi(x_n-1)
%}
% We can then update our matrix as follows
A_sub(1, 1:2) = [(2/3 + K*dx^2*x_sub(1)^2), -2/3]; % Update first row
A_sub(end, end-1:end) = [-2/3, (2/3 + K*dx^2*x_sub(end)^2)]; % Update last row

% Now can solve eigenvalue problem
[V, D] = eig(A_sub);
[d, ind] = sort(diag(D)); % sort eigenvalues
V_s = V(:, ind); % sort eigenvectors
eig_sols = d./(dx^2); % remember that these were for (eps*dx^2)*phi(x)
eig_sols = eig_sols(1:modes); % only want first five eigenvalues
V_s = V_s(:, 1:modes); % first five eigenvectors
V_sols = zeros(length(x_span), modes);

% Finally, we want to bootstrap eigenvectors
%{
  For this, we can use the previously derived forward and backward
  difference schemes. However, per the problem definition in this case we
  cannot assume that the dx term will go to zero. Thus, our value for
  phi(x_0) is as follows, noting that x_0^2 = (-L)^2 = L^2
    phi(x_0) = (4*phi(x_1) - phi(x_2))/(3 + 2*sqrt(L^2 - eps_0)*dx)

  Similarly, for the end point:
    phi(x_n) = (4*phi(x_n-1) - phi(x_n-2))/(3 + 2*sqrt(L^2 - eps_0)*dx)
%}
if plotting
    plt_labs = {}; % for storing plot labels
    f = figure("Name", "Direct Method");
    hold on; % to plot on same axis
end

for i=1:modes
    eig_bootstrap = [(4*V_s(1, i) - V_s(2, i))/ ...
                        (3 + 2*dx*sqrt(K*L^2 - eig_sols(i)));
                    V_s(:, i);
                    (4*V_s(end, i) - V_s(end-1, i))/...
                        (3 + 2*dx*sqrt(K*L^2 - eig_sols(i)))
                    ]; % create bootstrapped vector
    
    eig_bootstrap = abs(eig_bootstrap ./ ... % normalize bootstrapped vector
                    (sqrt(trapz(x_span, eig_bootstrap.*eig_bootstrap))));
    if plotting
        plt_labs{end + 1} = strcat("mode=", string(i)); % add label
        V_sols(:, i) = eig_bootstrap; % append results
        plot(x_span, V_sols(:, i)); % Plot results
    end
end

if plotting
    hold off;
    legend(plt_labs)
end

% ANSWERS
A7 = V_sols(:, 1); % First eigenvector
A8 = V_sols(:, 2); % Second eigenvector
A9 = V_sols(:, 3); % Third eigenvector
A10 = V_sols(:, 4); % Fourth eigenvector
A11 = V_sols(:, 5); % Fifth eigenvector
A12 = eig_sols; % eigenvectors

%% Problem 3 Adding a nonlinearity

modes = 2; % only considering two modes
eps_0 = 0.6;
L = 2; % shorter L
A0 = 0.3; % starting point for
gammas = [0.05 -0.05]; % gamma values
x_span = [-L:dx:L]; % new x_span
eigvec_nl = zeros(length(x_span), 4); % for storing solutions
eig_nl = []; % for storing eigenvalues

%{
  We can rewrite the function as a system of first order equations
    y_1 = phi, y_2 = y_1' = phi'
    ->
    y_1' = phi' = y_2
    y_2' = phi'' = (alpha*|y_1|^2 + K*x^2 - eps)*y_1
  
  We can implement a solution for this using a double shooting method
%}
if plotting
    f = figure("Name", "Double Shooting"); % create figure
    plt_labs = {}; % for storing labels
    hold on; % want all on the same plot
end

c = 1; % index for adding to eigenvector solutions
% BEGIN GAMMA LOOP
for g=1:length(gammas)
    gamma = gammas(g); % get value of gamma
    A_start = A0; % initial starting point for A
    eps_start = eps_0; % starting point for epsilon
    % BEGIN MODES LOOP
    for n=1:modes % number of modes to iterate through
        A = A0; % define A for loop
        dA = A_start / 100; % initial step for A shooting
        eps = eps_start;
        eps_step = eps_0 / 100;
        % BEGIN A SHOOTING LOOP
        for i=1:maxiter
            a_0 = [A A*sqrt(K*(L^2) - eps)]; % will be updated each loop with A
            [x, y] = ode45(@(x, y) nlShootingMethod(x, y, K, gamma, eps), ...
                                    x_span, a_0); % ode45
            % take want integral of vector squared - want this to equal 1 to satisfy integral -inf to inf phi^2 = 1
            sol_norm = trapz(x, y(:, 1).*y(:, 1));
            % if we are less than tolerance, we can move on to epsilon shooting
            if abs(sol_norm - 1) < tol
                break
            else
                A = A / sqrt(sol_norm);
            end
            % BEGIN EPSILON SHOOTING
            % initial condition for x using solved A
            x_0 = [A A*sqrt(K*(L^2) - eps)]; % Update initial guess
            [x_eps, y_eps] = ...
                ode45(@(x, y) nlShootingMethod(x, y, K, gamma, eps), ...
                                x_span, x_0); % ode45
            check = y_eps(end, 2) + sqrt(K*(L^2) - eps)*y_eps(end, 1); % boundary condition
            if abs(check) < tol
                break % End loop
            elseif (-1)^(n + 1) * check > 0 
                eps = eps + eps_step; % we are too low, need to be higher
            else % we are too high, need to decrease
                eps = eps - eps_step; % Update epsilon
                eps_step = eps_step / 2; % decrease eps_step
            end
        end
        eig_nl(end + 1) = eps; % add solution
        y_norm = abs(y_eps(:, 1) / sqrt(trapz(x_eps, y_eps(:, 1).*y_eps(:, 1)))); % normalize
        eigvec_nl(:, c) = y_norm; % append eigenvector
        c = c + 1; % update increment

        if plotting
            plot(x_eps, y_norm); % plot results
            plt_labs{end + 1} = strcat("mode=", string(n), ", gamma=", ...
                                        string(gamma));
        end
        eps_start = eps + 0.1; % increment epsilon
    end
end

if plotting
    legend(plt_labs); % add legend
    hold off;
end

% ANSWERS
A13 = eigvec_nl(:, 1); % first gamma = 0.5 vec
A14 = eigvec_nl(:, 2); % second gamma = 0.5 vec
A15 = eig_nl(1:2)'; % eigenvalues for gamma = 0.5
A16 = eigvec_nl(:, 3); % first gamma = -0.5 vec
A17 = eigvec_nl(:, 4); % second gamma = -0.5 vec
A18 = eig_nl(3:4)'; % eigenvalues for gamma = -0.5