clear;
close all;
clc;
plotting=true;

%% Default Parameters
nu = 0.001;
N = 256; % number of discretation points
L = 20; % size of discretation
m = 'FFT'; % will be solving exclusively using FFT method
x = linspace(-L/2, L/2, N+1);
x = x(1:N); % only consider first N points, periodic boundary conditions
y = x;
[x_grid, y_grid] = meshgrid(x, y); % Create grid structure for evaluating gaussians
t_span = 0:0.2:100; % time vector to use

% FFT specific parameters
kx = (2*pi/L)*[0:(N/2 - 1) (-N/2):-1]; % Define kx
kx(1) = 10^-6; % Update initial value to avoid divide by zero error
ky = kx; % We want to create a grid, so set kx = ky
[Kx, Ky] = meshgrid(kx, ky); % Create grid
K = Kx.^2 + Ky.^2; % Create the divisor for the FFT method

% Partial Matrices for solving systems
delta = abs(x(end) - x(end-1)); % delta value
A = generateLaplacian(N, delta); % laplacian matrix
B = generatePartialMatrices(N, delta, true); % partial x matrix
C = generatePartialMatrices(N, delta, false); % partial y matrix

% Gaussian function handle
%{
    Gaussian of the form
    A * exp(-(x-x0)
    Params:
    x - x coordinates to evaluate function at
    y - y coordinates to evaluate function at
    x_0 - x coordinate for centering gaussian
    y_0 - y coordinate for centering gaussian
    x_spread - x spread ( > 0)
    y_spread - y spread
    amp - amplitude adjustment (+ ccw, - cw for rotation)
%}
gaussian = @(x, y, x_0, y_0, x_spread, y_spread, amp) ...
        amp * exp(-((x - x_0).^2 .* (1 ./ (2 .* x_spread.^2)) + ...
                    (y - y_0).^2 .* (1 ./ (2 .* y_spread.^2))));

% Generate struct array with solve variables
solveVar.A = A; solveVar.B = B; solveVar.C = C;
solveVar.K = K; solveVar.N = N; solveVar.nu = nu;
solveVar.m = m;

%% Generate Solutions

for s=1:4 % iterate through three different scenarios
    
    if s == 1 % Oppositely Charged Vortices
        video_type = 'Opposite Charge'; % for video name
        w_grid = gaussian(x_grid, y_grid, 2, 0, 1, 3, 1) + ...
                gaussian(x_grid, y_grid, -2, 0, 1, 3, -1);        
    elseif s == 2 % Same charged vortices
        video_type = 'Same Charge'; % for video name
        w_grid = gaussian(x_grid, y_grid, 2, 0, 1, 3, 1) + ...
                gaussian(x_grid, y_grid, -2, 0, 1, 3, 1);
    elseif s == 3 % Pairs of charges
        video_type = 'Colliding Charges';
        % Below works, just want to mess with amplitude of second set
        w_grid = gaussian(x_grid, y_grid, 2, 2, 2, 0.5, -2) + ...
                gaussian(x_grid, y_grid, 2, -2, 2, 0.5, 2) + ...
                gaussian(x_grid, y_grid, -4, -2, 2, 0.5, 2) + ...
                gaussian(x_grid, y_grid, -4, 2, 2, 0.5, -2);
    else
        continue
    end
       
    % Solve
    w_0 = reshape(w_grid, N^2, []); % Reshape
    [t_res, y_res] = ode45(@(t, w) omegaTimestep(w, solveVar), t_span, w_0);
        
    if plotting
        f = figure("Name", strcat(m, ' AsCoolAsItGets - ', video_type)); % generate figure object
        grid off; % turn off grid
%         colormap("bone"); % set colormap
        
        for i=1:length(t_span)
            clf; % clear figure at each iteration
            y_i = reshape(y_res(i, :), N, []); % reshape to matrix form
            pcolor(x_grid, y_grid, y_i); % plot single image
            grid off; shading interp; % plot settings
            xlabel('x'); ylabel('y'); % Labels
            title(['Vorticity at time t=', num2str(t_span(i)), ' seconds']); % title
            movie_vector(i) = getframe(f); % get the frame
        end

        % Create the video writer and write out
        writer = VideoWriter(f.Name, 'MPEG-4'); % video writer name
        open(writer); % open the writer
        writeVideo(writer, movie_vector); % write video
        close(writer);

    end
end