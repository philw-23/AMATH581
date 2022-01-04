clear; close all; clc
plotting = true;

%% Define Initial Parameters

% Base Spectral Parameters
L = 20; % length of grid sides
N = 64; % discretation size
x = linspace(-L/2, L/2, N+1); % define x
x = x(1:N); % First N points for periodic BC
y = x; % same discretation in y
t_span = 0:0.5:4; % timespan
[x_grid, y_grid] = meshgrid(x, y); % grid structure

% Fourier domain parameters
kx = (2*pi/L)*[0:(N/2 - 1) (-N/2):-1]; % Define kx
kx(1) = 10^-6; % Update initial value to avoid divide by zero error
ky = kx; % We want to create a grid, so set kx = ky
[Kx, Ky] = meshgrid(kx, ky); % Create grid
K = -1*(Kx.^2 + Ky.^2); % Create the divisor for the FFT method

% Nonlinear Functions
lambda_A = @(U, V) 1 - (U.^2 + V.^2); % lambda function
omega_A = @(U, V, beta) -beta*(U.^2 + V.^2); % omega function

% Initial Condition Functions (Spirals)
u_init = @(m, x, y) tanh(sqrt(x.^2 + y.^2)).*cos(m*angle(x + 1i*y) - ...
                                            sqrt(x.^2 + y.^2));
v_init = @(m, x, y) tanh(sqrt(x.^2 + y.^2)).*sin(m*angle(x + 1i*y) - ...
                                            sqrt(x.^2 + y.^2));

% Define initial conditions
m = 1; % number of spirals
u = u_init(m, x_grid, y_grid);
v = v_init(m, x_grid, y_grid);
uf = fft2(u); % convert u to fourier space
vf = fft2(v); % convert v to fourier space
uvf_0 = [reshape(uf, N^2, []); reshape(vf, N^2, [])]; % concatenate initial condition

% Define General Parameters and build struct var
beta = 1; D1 = 0.1; D2 = 0.1;
svars.beta = beta; svars.D1 = D1; svars.D2 = D2; svars.N = N;
svars.lambda_A = lambda_A; svars.omega_A = omega_A; svars.lap = K;

%% Fourier Solutions

[tf_sol, yf_sol] = ode45(@(t, uvf) fftStep(t, uvf, svars), ...
                    t_span, uvf_0);

% SOLUTIONS
A1 = real(yf_sol); % Real portion of solution
A2 = imag(yf_sol); % Imaginary portion of solution

% Plot solution
if plotting
    plotSolutions(tf_sol, yf_sol, N, x_grid, y_grid, 'fft');
end

%% Define Cheb Parameters

N = 30; % cheb N
[Dmat, xc] = cheb(N); % calculate D matrix
D2mat = Dmat^2; % convert second derivative matrix
D2mat = (4/L^2) * D2mat; % scale to original discretation size
sCheb = size(D2mat, 2); % size of cheb matrix (<> N)
D2mat(1, :) = zeros(1, sCheb); % first row to zeros for BC
D2mat(end, :) = zeros(1, sCheb); % last row to zeros for BC
I = eye(length(D2mat)); % Identity matrix for kron
lapCheb = kron(D2mat, I) + kron(I, D2mat); % first term is dx^2, second dy^2

% Redefine grid
xc = L/2 .* xc; % scale to original discretation size
yc = xc; % same discretation in y
[x_grid, y_grid] = meshgrid(xc, yc); % create grid

% Redefine initial conditions
u = u_init(m, x_grid, y_grid);
v = v_init(m, x_grid, y_grid);
cheb_init = [reshape(u, sCheb^2, []); reshape(v, sCheb^2, [])]; % flatten to single vector

% Update struct var
svars.lap = lapCheb; svars.N = sCheb; % Note - size of cheb matrix <> N

%% Cheb Solutions

[t_cheb, y_cheb] = ode45(@(t, cheb) chebStep(t, cheb, svars), ...
                        t_span, cheb_init);

if plotting
    plotSolutions(t_cheb, y_cheb, sCheb, x_grid, y_grid, 'cheb');
end