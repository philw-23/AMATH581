clear;
close all;
clc;
plotting = true; % boolean for plotting

%% Exercise 1 - Generating Sparse Matrices

% Define variables
N = 8; % Number of steps
x_span = [-10:10]; % range of x values
dx = (max(x_span) - min(x_span)) / N; % step size
laplacian = generateLaplacian(N, dx);
derivativeX = generatePartialMatrices(N, dx, true);
derivativeY = generatePartialMatrices(N, dx, false);

% ANSWERS
A1 = full(laplacian);
A2 = full(derivativeX);
A3 = full(derivativeY);

%% Exercise 2 - FFT

% Load files
F = load('Fmat.mat').Fmat;
permvec = load('permvec.mat').permvec;

% Access Center Elements
N = size(F, 1); % size of square matrix
M = 20; % Block matrix (kernel) size
k = sqrt(length(permvec)); % side of square matrix representing number of MxM blocks 

% Access center elements
start_idx = (N - k * M) / 2 + 1; % index for start of center matrix
end_idx = start_idx + k * M - 1; % index for end of center matrix
Fcenter = F(start_idx:end_idx, start_idx:end_idx); % access center elements

% Get indice pairings for permutation
% ind2sub essentially gives us where each matrix block should live in the
% final result. These are in terms of a k x k matrix, but we can scale
% these values using M to access the full M x M block that we want to
% extract from Fcenter
[row_p, col_p] = ind2sub([k k], permvec);

% Get indice pairings for original Fcenter blocks
% This is similar to above, except in this case we want to place the
% extracted elements in order into our new center matrix. This will return
% k x k matrix indices in sequence 1:k^2, which is the desired order to
% place the extracted values in
[row_base, col_base] = ind2sub([k k], 1:length(permvec));

% Define permuted matrix and fill in values
Fpermuted = zeros(size(Fcenter));
for i=1:length(permvec) % iterate through permutation vector
    % Access desired elements from original matrix
    row_elems = (M * (row_p(i) - 1) + 1):(M * row_p(i)); % row elements to extract
    col_elems = (M * (col_p(i) - 1) + 1):(M * col_p(i)); % column elements to extract
    block_extract = Fcenter(row_elems, col_elems); % Extract values

    % Update permuted matrix
    row_update = (M * (row_base(i) - 1) + 1):(M * row_base(i)); % row elements to update
    col_update = (M * (col_base(i) - 1) + 1):(M * col_base(i)); % col elements to update
    Fpermuted(row_update, col_update) = block_extract; % Make updates
end

% Recenter matrix
Ftransformed = F;
Ftransformed(start_idx:end_idx, start_idx:end_idx) = Fpermuted;

% Unshift and transform
Fshifted = ifftshift(Ftransformed); % shift fft
f_no_fft = ifft2(Fshifted); % inverse transform on fft

% ANSWERS
A4 = abs(Ftransformed); % Decrypted fourier matrix
A5 = abs(f_no_fft); % recovered image matrix

if plotting % Plot resulting image
    figure('Name', 'Decrypted Image');
    set(gcf, 'colormap', gray);
    imagesc(A5);
end