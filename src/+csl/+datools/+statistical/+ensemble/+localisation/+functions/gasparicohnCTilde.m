function Ctilde = gasparicohnCTilde(n, r, d, t, y, m, k, theta)

rsfun = @(r) repmat(r, n/numel(r), 1);

r = rsfun(r);

if nargin < 8
    theta = 1;
end

CtildeD = zeros(n, 1);

if isscalar(r)
    r = r*ones(n,1);
end

for l = 1:n
    if ~isempty(m)
        % NEW HOTNESS
        ki = d(t, y, k, l)/r(k);
        kj = d(t, y, k, l)/r(l);
        CtildeD(l) = m(gc(ki, theta), gc(kj, theta));
    else
        ki = d(t, y, k, l)/r(k);
        CtildeD(l) = gc(ki, theta);
    end
end

Ctilde = spdiags(CtildeD, 0, n, n);

end

function g = gc(k, theta)

if k <= theta
    g = 1 - (5/3)*k^2 + (5/8) * k^3 + (1/2)*k^4 - (1/4)*k^5;
elseif k <= 2*theta
    g = 4 - 5*k + (5/3)*k^2 + (5/8)*k^3 - (1/2)*k^4 + (1/12)*k^5 - (2/(3*k));
else
    g = 0;
end

end
