%{
    A function for generating a laplacian matrix
    Inputs:
        N - integer representing the size of the matrix discretation
        dx - step size for discretation
    Returns:
        L - N^2 x N^2 matrix representing the laplacian conversion matrix
%}
function L = generateLaplacian(N, dx)

    % Generate items for blocks
    % Vectors are generated as length N
    base_ones = ones(N, 1); % base ones vector to use
    d_main = -4 * base_ones; % Main diagonal of block
    d_off_one = base_ones; % initialize off diagonal
    d_off_one(end) = 0; % pattern is N-1 ones, then one 0
    d_off_nminus = base_ones; % n-1 off diagonal
    d_off_nminus(2:end) = 0; % pattern is one 1, N-1 zeros

    % Generate Matrix for input to spdiags
    bin_mat = repmat([base_ones, ... % -1*(N^2 - N) diag
                    base_ones, ... % -1*N diag
                    d_off_nminus, ... % -1 * (N - 1) diag
                    d_off_one, ... % -1 diag
                    d_main, ... % main diagonal
                    d_off_one, ... % +1 diag
                    d_off_nminus, ... % N - 1 diag
                    base_ones, ... % N diag
                    base_ones % N^2 - 1 diag
                ], N, 1); % repeat N times

  
    % Generate vector of diagonal locations
    diags = [-(N^2 - N) -N -(N - 1) -1 0 1 (N - 1) N (N^2 - N)];

    %{
        Note: for diagonals < 0 when num_rows >= num_cols we have no needed 
        action as the top values will be selected first. However, for
        diagonals > 0 the bottom values will be selected first so we need
        to shift the bin values down based on the diagonal. This can be
        accomplished using the circshift command
    %}
    for i=1:length(diags)
        d = diags(i); % Helper variable
        if d > 0 % if we are upper diagonal we wan to slide values
            bin_mat(:, i) = circshift(bin_mat(:, i), d);
        end
    end

    % Generate matrix
    L = spdiags(bin_mat, diags, N^2, N^2);

    % Divide by dx^2 per second order central difference formula
    L = L ./ (dx^2);

end