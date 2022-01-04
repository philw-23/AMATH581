clear;
close all;
clc;
plotting=true;

%% Initialize Parameters

% Default model parameters
N = 16; % spacial discretization
L = 2*pi; % going from -pi to pi
x = linspace(-L/2, L/2, N+1); % discretize with N+1 points
x = x(1:N); % take first N points
y = x; z = x; % same discretation in other dimensions
[x_grid, y_grid, z_grid] = meshgrid(x, y, z); % create grid
t_span = 0:0.5:4; % timespan of interest
A = [-1 -1 -1]; % Define A parameters
B = -1 .* A; % Define B parameters

% Non-Linear Variables
% Sin terms are static and can be calculated outside of function
sin_terms = (A(1) .* sin(x_grid).^2 + B(1)) .* ...
            (A(2) .* sin(y_grid).^2 + B(2)) .* ...
            (A(3) .* sin(z_grid).^2 + B(3));

% FFT conditions
% Values and variables needed for calculations in fourier space
% Applies to linear potion of problem: 1/2 * laplacian(phi)
%{
    We can derive the laplacaian in terms of fourier vectors Kx, Ky, and Kz
    using the derivative property of the FFT
    
    laplacian = d^2/dx^2 + d^2/dy^2 + d^2/dz^2
    Fourier space derivative property: f^(n) = (ik)^n
    therefore laplacian = (iKx)^2 + (iKy)^2 + (iKz)^2
        = -Kx^2 - Ky^2 - Kz^2
        = -(Kx^2 + Ky^2 + Kz^2)
%}
kx = (2*pi/L)*[0:(N/2 - 1) (-N/2):-1]; % Define kx
kx(1) = 10^-6; % Update initial value to avoid divide by zero error
ky = kx; kz = kx; % same discretation in other dimensions 
[Kx, Ky, Kz] = meshgrid(kx, ky, kz); % Create grid
K = -(Kx.^2 + Ky.^2 + Kz.^2); % Laplacian operator in fourier space

% Define struct vector with variables needed to solve
vars.N = N; vars.L = L; vars.x_grid = x_grid; vars.y_grid = y_grid; 
vars.z_grid = z_grid; vars.K = K; vars.sin_terms = sin_terms;

%% Part 1 - Cosine Initial Conditions

% Define initial condition
phi_0_grid = cos(x_grid).*cos(y_grid).*cos(z_grid);
phi_0_grid_fft = fftn(phi_0_grid); % Take FFT
phi_0_fft = reshape(phi_0_grid_fft, N^3, []); % flatten for ode45

% solve using ode45
[t, y_cosine] = ode45(@(t, phi) bec3DTimestep(t, phi, vars), ...
                        t_span, phi_0_fft);

% SOLUTIONS
A1 = real(y_cosine); % real portion of solution
A2 = imag(y_cosine); % imaginary portion of solution

if plotting
    % Generate gifs
    filename = 'BEC - Cosine Initial Condition'; % create filename
    plotSolutions(t, y_cosine, x_grid, y_grid, z_grid, N, 'gif', filename);
    % Generate time progressions
    plotSolutions(t, y_cosine, x_grid, y_grid, z_grid, N, 'prog', filename);
    plotSolutions(t, y_cosine, x_grid, y_grid, z_grid, N, 'contour', filename);
end

%% Part 2 - Sine Initial Conditions

% Define initial condition
phi_0_grid = sin(x_grid).*sin(y_grid).*sin(z_grid);
phi_0_grid_fft = fftn(phi_0_grid); % Take FFT
phi_0_fft = reshape(phi_0_grid_fft, N^3, []); % flatten for ode45

% solve using ode45
[t, y_sine] = ode45(@(t, phi) bec3DTimestep(t, phi, vars), ...
                    t_span, phi_0_fft);

% SOLUTIONS
A3 = real(y_sine); % real portion of solution
A4 = imag(y_sine); % imaginary portion of solution

if plotting
    % Generate Gifs
    filename = 'BEC - Sine Initial Condition'; % create filename
    plotSolutions(t, y_sine, x_grid, y_grid, z_grid, N, 'gif', filename);
    % Generate time progressions
    plotSolutions(t, y_sine, x_grid, y_grid, z_grid, N, 'prog', filename);
    plotSolutions(t, y_sine, x_grid, y_grid, z_grid, N, 'contour', filename);
end

% Close all figures
close all; % close figures