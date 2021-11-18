clear;
close all;
clc;
plotting=false; % boolean for plotting

%% Vorticity-Streamfunction Parameters

%{
    The purpose of this exercise is to add timestepping to the
    streamfunctions that we have previously been working with and testing
    out different methods. We will be discretizing the following two
    streamfunction equations

    w_t + [phi, w] = nu * laplacian(w)
    w/ [phi, w] = phi_x * w_y - phi_y * w_x

    laplacian(phi) = w
%}
% Define initial parameters
nu = 0.001; % nu parameter
N = 64; % number of discretations
L = 20; % length of x, y

% Generate discretation for normal space
x = linspace(-L/2, L/2, N + 1); % x vector
y = linspace(-L/2, L/2, N + 1); % y vector
% NOTE: We are assuming periodic boundary conditions, thus we generate x
% and y vectors with N + 1 points and drop the last point.
x = x(1:N); % Take first 64 points
y = y(1:N); % Take first 64 points
[x_grid, y_grid] = meshgrid(x, y); % generate discretized grid of x, y
t_span = 0:0.5:4; % timespan for problem

%{
    We also will need to solve equation 2 in fourier space, so we will
    perform the necessary discretation now to pass to the solver. To gain
    back the solution in normal space, we can take the inverse fourier
    transform of the result

    laplacian(phi) = w ==> phi = -w / kx.^2 + ky.^2
%}
kx = (2*pi/L)*[0:(N/2 - 1) (-N/2):-1]; % Define kx
kx(1) = 10^-6; % Update initial value to avoid divide by zero error
ky = kx; % We want to create a grid, so set kx = ky
[Kx, Ky] = meshgrid(kx, ky); % Create grid
K = Kx.^2 + Ky.^2; % Create the divisor for the FFT method

%{
    We can rewrite the system as follows:
        w_t = nu * laplacian(w) - [phi, w]
        w_t = nu * laplacian(w) - phi_x * w_y + phi_y * w_x

    From our previous matrix discretations of phi, we know that this can be
    written as follows where A, B, and C are discretized partial matrices.
    A represents the laplacian operator, B represents the partial with
    respect to x and C represents the partial with respect to y
        w_t = nu * A*w - (B*phi)*(C*w) + (C*phi)*(B*w)
%}
% Define our partial matrices
delta = abs(x(end) - x(end-1)); % delta value
A = generateLaplacian(N, delta); % laplacian matrix
B = generatePartialMatrices(N, delta, true); % partial x matrix
C = generatePartialMatrices(N, delta, false); % partial y matrix
% NOTE: to avoid singularity in the laplacian matrix, we must set the first
% value to be <> -4 when using LU decomp and backslash to solve for phi. 
% A can stay the same for solving the full equation
AMod = A;
AMod(1, 1) = 2 ./ (delta^2); % Set to two per problem definition, but other values could be used
[Ld, U, P] = lu(AMod); % LU Decomp

%{
    Finally, want to define our initial time condition, evaluate on the 
    grid and reshape to be a long vector for the vorticity function w. We
    can then use this to solve for phi to use in our ode solvers using the
    following equation

        laplacian(phi) = A*phi = w
%}
gaussian_bump = @(x, y) exp(-x.^2 - (y.^2 ./ 20));
w_grid = gaussian_bump(x_grid, y_grid);
w_0 = reshape(w_grid, N^2, 1);

%% Solve Using Various Methods

% Will test multiple methods to solve for phi
methods = {'Backslash', 'LU Decomp', 'bicgstab', 'gmres', 'FFT'};
solutions = zeros(9, 4096, length(methods)); % For storing solution vectors
times = zeros(length(methods), 1); % Storage for timing matrices
maxiterB = 1000; % max iterations for bicgstab
restart = 10; % restart iterations for gmres
maxiterG = maxiterB / restart; % iteration sets for gmres
tol = 1e-8; % tolerance for bicgstab and gmres

% Create a struct object for items we need for different solvers
solveVar.A = A; solveVar.B = B; solveVar.C = C; solveVar.AMod = AMod;
solveVar.K = K; solveVar.N = N; solveVar.nu = nu; 
solveVar.maxiterB = maxiterB; solveVar.maxiterG = maxiterG; 
solveVar.tol = tol; solveVar.restart = restart;
solveVar.Ld = Ld; solveVar.U = U; solveVar.P = P;

for i=1:length(methods)
    
    % SOLVE FOR PHI USING MULTIPLE METHODS
    m = methods{i};
    disp(m);
    solveVar.m = m;
    tic; % start time for execution
    
    % SOLVE ODE
    [t_res, y_res] = ode45(@(t, w) omegaTimestep(w, solveVar), ...
                       t_span, w_0); 

    % APPEND RESULTS
    times(i) = toc; % append solve times
    solutions(:, :, i) = y_res; % Add results
    
    % Generate gifs for reference if we are plotting
    if plotting

        f = figure("Name", m); % generate figure object
        grid off; % turn off grid
        filename = strcat(m, '.gif'); % generate gif filename

        for j=1:length(t_span) % iterate through all times
            y_i = reshape(y_res(j, :), N, []); % reshape to matrix form
            pcolor(y_i); % plot single image
            shading interp; % remove grid blocks
            xlabel('x'); ylabel('y'); % axis labels
            title(['Vorticity at time t = ', ...
                num2str(t_span(j)), ' seconds']); % title

            % Get frame and convert to image
            frame = getframe(f); % get the frame for f
            image = frame2im(frame); % convert frame to image
            [imind, cmap] = rgb2ind(image, 256); % convert to indexed image

            % Write out to gif file
            if j == 1 % generate initial image write
                imwrite(imind, cmap, filename, 'gif', 'LoopCount', inf);
            else % add additional frames
                imwrite(imind, cmap, filename, 'gif', 'WriteMode', 'append');
            end
        end
    end
end

% ANSWERS
A1 = solutions(:, :, 1); % Backslash
A2 = solutions(:, :, 2); % LU Decomp
A3 = solutions(:, :, 3); % bicgstab
A4 = solutions(:, :, 4); % gmres
A5 = solutions(:, :, 5); % FFT