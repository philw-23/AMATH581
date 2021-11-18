% RHS ode function for the votrticity defined as:
% w_t = nu * A*w - (B*phi)*(C*w) + (C*phi)*(B*w)
% Note that many of these operators will be matrices so matrix operations
% need to be considered
function dwdt = omegaTimestep(w, vars)
    
    if strcmp(vars.m, 'Backslash')
        % Solve for phi_0 using the backslash method
        phi = vars.AMod \ w; % solving laplacian(phi) = w
    elseif strcmp(vars.m, 'LU Decomp')
        % Solve for phi_0 using LU decomposition
        y = vars.Ld \ (vars.P^(-1) * w); % solve for y 
        phi = vars.U \ y; % solve for phi
    elseif strcmp(vars.m, 'bicgstab')
        % solve for phi_0 using bicgstab method
        % bi-conjugate stabilized gradient
        [phi, ~, ~, ~] = bicgstab(vars.AMod, w, ...
                            vars.tol, vars.maxiterB);
    elseif strcmp(vars.m, 'gmres')
        % solve for phi_0 using gmres method
        % generalized minimum residual method
        [phi, ~, ~, ~] = gmres(vars.AMod, w, vars.restart, ...
                            vars.tol, vars.maxiterG);
    elseif strcmp(vars.m, 'FFT') 
        w_grid = reshape(w, vars.N, []); % Reshape to grid form
        w_fft = fft2(w_grid); % Take FFT of initial condition
        phi = real(ifft2(-w_fft./vars.K)); % Compute phi in fft space
        phi = reshape(phi, vars.N^2, []); % Reshape to vector form
    end

    dwdt = vars.nu.*(vars.A*w) - ...
            (vars.B*phi).*(vars.C*w) + ...
            (vars.C*phi).*(vars.B*w);

end

