function [xz, xzPlot] = sineWave(N, xLim, scale, sigma)
if nargin<2
    xLim = [0, 3*pi];
if nargin<3
    scale = 0.2;
end
if nargin<4
    sigma = 0.25;
end
xLeft = xLim(1);
xRight = xLim(end);
x = linspace(xLeft, xRight, N);

map = @(x,sig,n) exp((-scale*x) .* sin(x)) + sig*randn(1,n);
z = map(x,sigma, N);

xz = [x;z];
xzPlot = xz;
end
