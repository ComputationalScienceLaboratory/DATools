function rhoHt = gc(t, y, r, d, H, m)

if nargin < 6
   m = @(ri, rj) (ri + rj)/2; 
end

n = size(H, 2);

I2 = 1:n;

I1 = find(sum(abs(H), 1)).';

rhoHt = zeros(n, numel(I1));
if numel(r) ~= n
    r = repmat(r, n/numel(r), 1);
end
r = r.';

for jr = 1:numel(I1)
    
    j = I1(jr);
    
    ks1 = d(t, y, I2, j)./r(j);
    ks2 = d(t, y, I2, j)./r(I2);
    
    
    rhoHt(:, jr) = m(gcf(ks1).', gcf(ks2).');
end

end

function g = gcf(k)

theta = 1;

g = zeros(size(k));

mask1 = k <= theta;
mask2 = and((k > theta), (k <= 2*theta));

km1 = k(mask1);
km2 = k(mask2);

g(mask1) = ones(size(km1)) - (5/3)*(km1.^2) + (5/8) * (km1.^3) + (1/2)*(km1.^4) - (1/4)*(km1.^5);
g(mask2) = 4*ones(size(km2)) - 5*km2 + (5/3)*(km2.^2) + (5/8)*(km2.^3) - (1/2)*(km2.^4) + (1/12)*(km2.^5) - (2./(3*km2));

end