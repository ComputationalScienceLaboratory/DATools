function Ctilde = cutoffCTilde(t, y, H, r, d, k)

%rsfun = @(r) repmat(r, n/numel(r), 1);

%r = rsfun(r);

ki = d(t, y, k, H) / r;
CtildeD = cut(ki).';

Ctilde = spdiags(CtildeD, 0, numel(H), numel(H));

end

function g = cut(k)


g = zeros(size(k));

mask = k <= sqrt(pi/2);

%km1 = 1 - k(mask)/2;

km1 = ones(size(k(mask)));

g(mask) = ones(size(km1));

end
