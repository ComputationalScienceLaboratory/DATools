function y = rosEu1step(f, tspan, y0, J)

n = length(y0);
h = tspan(2);

feval = f(tspan(1), y0);

Jeval = eye(n) - h*J(tspan(1), y0);
% Jeval = eye(n);

y = y0 + h*(Jeval\feval);

end