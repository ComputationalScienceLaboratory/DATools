function [xz, xzPlot] = normalGamma(N, a, b, sigma)
if nargin<2
    a = 2;
if nargin<3
    b = 2;
end
if nargin<4
    sigma = 0.25;
end

x = gamrnd(a, b, [1, N]);
map = @(x, sig, n) x + sig*randn(1, n);

z = map(x, sigma, N);

xz = [x; z];
xzPlot = xz;
end