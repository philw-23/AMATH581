%{
    A function for generating first order dt^2 partialx and partialy
    matrices

    Inputs:
        N - integer representing the size of the matrix discretation
        delta - step size for discretation
        partialx - boolean indicating x or y partial
    Returns:
        P - N^2 x N^2 matrix representing the partial conversion matrix
%}
function P = generatePartialMatrices(N, delta, partialx)

    % Generate items for blocks
    % Vectors are generated as length N
    base_ones = ones(N, 1); % base ones vector to use
    d_neg = -1 * base_ones; % negative ones vector
  
    % Generate vector of diagonal locations
    if partialx % case for x derivative
        bin_mat = repmat([base_ones, ...
                        d_neg, ...
                        base_ones, ...
                        d_neg
                    ], N, 1); % repeat N times
        diags = [-(N - 1)*N -N N (N - 1)*N];
    else
        base_ones_y = base_ones;
        base_ones_y(2:end) = 0; % Pattern is one 1, N-1 zeros
        base_neg_ones = d_neg;
        base_neg_ones(end) = 0; % pattern is N-1 -1's, 1 zero

        bin_mat = repmat([base_ones_y, ...
                        base_neg_ones, ...
                        -1 * base_neg_ones, ...
                        -1 * base_ones_y
                    ], N, 1); % repeat N times

        diags = [-(N - 1) -1 1 (N - 1)];
    end

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
    P = spdiags(bin_mat, diags, N^2, N^2);

    % Divide by 2*delta per first order central difference formula
    P = P ./ (2*delta);

end