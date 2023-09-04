function rhoHt = gauss(t, y, r, d, H, m)

if nargin < 6
    m = @(ri, rj) (ri + rj) / 2;
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

    ks1 = d(t, y, I2, j) ./ r(j);
    ks2 = d(t, y, I2, j) ./ r(I2);


    rhoHt(:, jr) = m(exp(-(1 / 2)*(ks1.^2)).', exp(-(1 / 2)*(ks2.^2)).');
end

end