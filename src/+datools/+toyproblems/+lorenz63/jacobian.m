% this is the jacobian of the rhs of the lorenz63 equations
%
function J = jacobian(~, y, sigma, rho, beta)
if nargin<3
    sigma = 10;
    rho = 28;
    beta = 8/3;
end
J = [-sigma  ,  sigma,  0   ; ...
    rho-y(3), -1    , -y(1); ...
    y(2)    ,  y(1) , -beta];

end