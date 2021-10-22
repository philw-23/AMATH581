function dYdx = oscillatorShoot(x, y, K, eps)
    % Note: x is essentially t in this
    dYdx = [y(2); (K*(x^2) - eps)*y(1)];

end