% Nonlinear shooting function derivative
function dYdt = nlShootingMethod(x, y, K, gamma, eps) 
    dYdt = [y(2);
            (gamma*abs(y(1))^2 + K*x^2 - eps)*y(1)];
end