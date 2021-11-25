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
t_span = 0:0.2:100; % time vector to use for solving/plotting

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

for s=1:6 % iterate through six different scenarios
   
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
        % Swapping the rotation/amplitude of the left two points causes the
        % particle pairs to charge towards each other and then drift apart.
        w_grid = gaussian(x_grid, y_grid, 2, 2, 2, 0.5, -2) + ...
                gaussian(x_grid, y_grid, 2, -2, 2, 0.5, 2) + ...
                gaussian(x_grid, y_grid, -4, -2, 1.5, 0.5, 2) + ...
                gaussian(x_grid, y_grid, -4, 2, 1.5, 0.5, -2);

    elseif s == 4 % Random set charges
        video_type = 'Random Charges';
        r_seed = rng(10); % set seed
        num_charge = 15; % 15 charges
        center_loc = L/2 - 2; % Maximum x and y center
        x_center = -center_loc + 2 * center_loc * rand([1 num_charge]); % x centers
        y_center = -center_loc + 2 * center_loc * rand([1 num_charge]); % y centers
        amp = -4 + (4+4)*rand([1 num_charge]); % amplitude between -4 and 4
        x_spreads = 0.5 + (2 - 0.5) * rand([1 num_charge]); % x spreads
        y_spreads = 0.5 + (2 - 0.5) * rand([1 num_charge]); % x spreads
        w_grid = zeros(size(x_grid)); % for adding random gaussians to

        for i=1:num_charge % iterate through number of charges
            w_grid = w_grid + gaussian(x_grid, y_grid, x_center(i), ...
                                y_center(i), x_spreads(i), y_spreads(i), ...
                                amp(i)); % add random gaussian to grid
        end

    elseif s >= 5 % Custom Radial Charges 
        if s == 5
            video_type = 'Radial Charges (Constant Amp)';
            amp_f = @(r_0, r) 2 * (-1)^r; % amplitude function
        elseif s == 6
            video_type = 'Radial Charges (Variable Amp)';
            amp_f = @(r_0, r) r_0 * (-1)^r; % amplitude function
        end


        r_max = 8; % Max radius to go to
        w_grid = zeros(size(x_grid)); % for adding radial gaussians
        
        x_unit = [cosd(0) cosd(30) cosd(60) cosd(90) cosd(120) ...
                        cosd(150) cosd(180) cosd(210) cosd(240) ...
                        cosd(270) cosd(300) cosd(330)]; % x unit circle locs
        y_unit = [sind(0) sind(30) sind(60) sind(90) sind(120) ...
                        sind(150) sind(180) sind(210) sind(240) ...
                        sind(270) sind(300) sind(330)]; % y unit circle locs
        x_spread = 0.5; % static x_spread
        y_spread = 0.5; % static y_spread

        for r=1:r_max % iterate through radiuses
            x_centers = r .* x_unit; % x centers
            y_centers = r .* y_unit; % y centers
            amp = amp_f(r, r); % amplitude

            for i=1:length(x_centers) % iterate through charges
                w_grid = w_grid + gaussian(x_grid, y_grid, x_centers(i), ...
                                    y_centers(i), x_spread, y_spread, ...
                                    amp); % add in charges
            end
        end
    end
       
    % Solve
    w_0 = reshape(w_grid, N^2, []); % Reshape
    [t_res, y_res] = ode45(@(t, w) omegaTimestep(w, solveVar), t_span, w_0);
        
    if plotting
        f = figure("Name", strcat([m, ' AsCoolAsItGets - ' , video_type])); % generate figure object
        grid off; colormap("turbo"); % set colormap, turn off grid
        hold on; % Don't create new figure in for loop

        for i=1:length(t_span)
            y_i = reshape(y_res(i, :), N, []); % reshape to matrix form
            pcolor(x_grid, y_grid, y_i); % plot single image
            grid off; shading interp % plot settings
            xlabel('x'); ylabel('y'); % Labels
            title(['Vorticity at time t=', num2str(t_span(i)), ' seconds']); % title
            movie_vector(i) = getframe(f); % get the frame
            clf(f); % Clear figure
        end

        % Create the video writer and write out
        writer = VideoWriter(f.Name); % video writer name
        writer.FrameRate = 10; % lower frame rate to 10
        open(writer); % open the writer
        writeVideo(writer, movie_vector); % write video
        close(writer);
        close(f); % close figure
        clear f; % clear figure

    end
end