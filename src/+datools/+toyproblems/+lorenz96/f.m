% this is the rhs (f) of Lorenz 96 dynamical system
%
function dy = f(~, y, N, F)
if nargin<3
    N = 40; % default
    F = 8; % creates chaotic map
end
if nargin<4
    F = 8; % creates chaotic map
end

dy = zeros(N, 1);
for i = 1:N
    dy(i) = (y(mod(i, N) + 1) - y(mod(i - 2 - 1, N) + 1)) * y(mod(i - 1 - 1, N) + 1) - y(i) + F;
end

end