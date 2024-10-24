% this is the rhs of Lorenz 63 dynamical system
%
function dy = f(~, y, sigma, rho, beta)
if nargin<3
    sigma = 10;
    rho = 28;
    beta = 8/3;
end
dy = [sigma * (y(2, :) - y(1, :)); ...
    rho * y(1, :) - y(2, :) - y(1, :) .* y(3, :); ...
    y(1, :) .* y(2, :) - beta * y(3, :)];

end