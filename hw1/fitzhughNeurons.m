% Defines a function to use to solve ode problem for coupled neurons
% inputs: 
%       t = time vector (has to exist, but not always passed)
%       y = y vector
%       a, b, c, d, I = parameters governing interaction between equations
%       (**Note: a, d, and y are vector inputs**)
% returns: 
%       dYdt = column vector of four differential equations to pass to ode solver:
%               [dv1; dv2; dw1; dw2]
function dYdt = fitzhughNeurons(t, y, a, b, c, d, I)

    % Define function pieces based on order we pass them in
    v1 = y(1);
    v2 = y(2);
    w1 = y(3);
    w2 = y(4);

    % Define derivatives
    % NOTE: a, d are vectors
    dv1 = -v1^3 + (1 + a(1))*v1^2 - a(1)*v1 - w1 + I + d(1)*v2;
    dv2 = -v2^3 + (1 + a(2))*v2^2 - a(2)*v2 - w2 + I + d(2)*v1;
    dw1 = b*v1 - c*w1;
    dw2 = b*v2 - c*w2;

    % Concatenate the results
    dYdt = [dv1; dv2; dw1; dw2];

end