function cheb_step = chebStep(t, cheb_val, vars)
    
    % Split U and V
    % Note - due to matrix sizes we don't reshape
    u = cheb_val(1:vars.N^2); % U is first
    v = cheb_val((vars.N^2 + 1):end); % V second

    % Calculate lambda and omega
    lambda_A = vars.lambda_A(u, v);
    omega_A = vars.omega_A(u, v, vars.beta);

    % Compute steps
    u_nl_term = lambda_A .* u  - omega_A .* v;
    v_nl_term = omega_A .* u + lambda_A .* v;
    u_step = u_nl_term + (vars.D1 .* vars.lap * u);
    v_step = v_nl_term + (vars.D2 .* vars.lap * v);

    % Flatten and return
    cheb_step = [u_step; v_step];

end