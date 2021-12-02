function uvf_step = fftStep(t, uvf, vars)

    % Split u and v
    uf = uvf(1:vars.N^2); % U is first N^2 points
    vf = uvf((vars.N^2 + 1):end); % V is next N^2 points
    
    % Reshape U and V
    uf_grid = reshape(uf, vars.N, []); % reshape u to grid form
    vf_grid = reshape(vf, vars.N, []); % reshape v to grid form

    % Convert to real space
    u_grid = real(ifft2(uf_grid)); % convert out of fourier space
    v_grid = real(ifft2(vf_grid)); % convert v out of fourier space
    
    % Evaluate non-linearities
    lambda = vars.lambda_A(u_grid, v_grid); 
    omega = vars.omega_A(u_grid, v_grid, vars.beta);
    u_nl = lambda .* u_grid - omega .* v_grid; % non-linear portion for U_t
    v_nl = omega .* u_grid + lambda .* v_grid; % non-linear portion for V_t
    u_nl_fft = fft2(u_nl); % convert u non-linear to fourier space
    v_nl_fft = fft2(v_nl); % convert v non-linear to fourier space

    % Evaluate linear operators in fourier space
    u_l_fft = vars.D1 .* vars.lap .* uf_grid; % laplacian U in fourier
    v_l_fft = vars.D2 .* vars.lap .* vf_grid; % laplacian V in fourier

    % Combine fourier results for next step
    uf_next = reshape(u_nl_fft + u_l_fft, vars.N^2, []); % combine and reshape U
    vf_next = reshape(v_nl_fft + v_l_fft, vars.N^2, []); % combine and reshape V
    uvf_step = [uf_next; vf_next];

end