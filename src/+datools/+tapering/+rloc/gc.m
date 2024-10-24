function Ctilde = gc(y, r, d, H, k)

%rsfun = @(r) repmat(r, n/numel(r), 1);

%r = rsfun(r);

ki = d(y, k, H) / r;
CtildeD = gcf(ki).';

%Ctilde = spdiags(CtildeD, 0, numel(H), numel(H));

Ctilde = diag(CtildeD);

end

function g = gcf(k)

% This is the solution to the equation
%
%   \int_0^\infty gc(x, \theta) dx = \int_0^\infty e^{- (x^2) / 2} dx
%
% and is roughly 1.77883918308445544346236005321
theta = (3 * sqrt(2*pi)) / (7 - 4 * log(2));

% Another option would be to find the argmin of
%
%   \int_0^(2\theta) [gc(x, \theta) - e^{- (x^2) / 2}]^2 dx
%
% which is roughly 1.76479800596413372026916022151
% If we instead integrate to infinity, the answer is roughly
% 1.76481571580742002680608493392
%
% We will use the first one as it has a closed form solution and the
% numbers are roughly the same, so it should not make that much of a
% difference.

g = zeros(size(k));

mask1 = k <= theta;
mask2 = and((k > theta), (k <= 2 * theta));

km1 = k(mask1) / theta;
km2 = k(mask2) / theta;

g(mask1) = ones(size(km1)) - (5 / 3) * (km1.^2) + (5 / 8) * (km1.^3) + (1 / 2) * (km1.^4) - (1 / 4) * (km1.^5);
g(mask2) = 4 * ones(size(km2)) - 5 * km2 + (5 / 3) * (km2.^2) + (5 / 8) * (km2.^3) - (1 / 2) * (km2.^4) + (1 / 12) * (km2.^5) - (2 ./ (3 * km2));

end
