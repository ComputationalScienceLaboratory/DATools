function rhoHt = cutoff(t, y, r, d, H, m)

if nargin < 6
   m = @(ri, rj) (ri + rj)/2; 
end

n = size(H, 2);

I2 = 1:n;

I1 = find(sum(abs(H).', 2));

%I1 = 1:size(H, 1);

%rhoHt = zeros(n, numel(I1));

rhoHt = zeros(n, size(H.', 2));
if numel(r) ~= n
    r = repmat(r, n/numel(r), 1);
end
r = r.';

for jr = 1:numel(I1)
    
    j = I1(jr);
    
    ks1 = d(t, y, I2, j)./r(j);
    ks2 = d(t, y, I2, j)./r(I2);
    
    
    %rhoHt(:, jr) = m(gcf(ks1).', gcf(ks2).')*(H(:, j).');
    
    rhoHt = rhoHt + m(cut(ks1).', cut(ks2).')*(H(:, j).');
end

end

function g = cut(k)


g = zeros(size(k));

mask = k <= 2;

%km1 = 1 - k(mask)/2;

km1 = ones(size(k(mask)));

g(mask) = ones(size(km1));

end
