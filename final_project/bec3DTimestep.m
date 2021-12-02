%{
Function to solve a timestep of the 3D Bose Einstein Condensation.
    i*phi_t + (1/2)*laplacian(phi) - |phi|^2*phi + ...
    (A1*sin^2(x) + B1)(A2*sin^2(y) + B2)(A3*sin^2(z) + B3)*phi = 0

Solving for phi_t, this equates to:
    phi_t = (1/i)*(|phi|^2 * phi - (1/2)*laplacian(phi) - ...
    (A1*sin^2(x) + B1)(A2*sin^2(y) + B2)(A3*sin^2(z) + B3)*phi)

Inputs:
    t - timespan of consideration
    phi - condition at the beginning of each timestep
    vars - struct array containing variables needed for solving equations
%}
function phi_step = bec3DTimestep(t, phi, vars)

    % Reshape input to grid form
    phi_res_fft = reshape(phi, vars.N, vars.N, []); % reshape to grid form
    
    % Linear portion
    lin_terms_fft = -0.5 .* vars.K .* phi_res_fft; % laplacian operator

    % Non linear portions
    phi_real_grid = ifftn(phi_res_fft); % take inverse FFT
    nl_sin_real = -phi_real_grid .* vars.sin_terms; % multiply by phi for sin terms
    nl_sin_fft = fftn(nl_sin_real); % Take FFT of sin terms
    phi_nl_real = phi_real_grid .* abs(phi_real_grid).^2; % Calculate non linear phi terms
    phi_nl_fft = fftn(phi_nl_real); % Take FFT of non-linear terms

    % Combine results
    phi_step_grid = (nl_sin_fft + phi_nl_fft + lin_terms_fft) ./ 1i; % Combine results
    phi_step = reshape(phi_step_grid, vars.N^3, []); % reshape to flat form

end