function rhoHt = gasparicohnrho_small(n, r, d, t, y, m, theta, H)

rsfun = @(r) repmat(r, n/numel(r), 1);

r = rsfun(r);

if nargin < 7
    theta = 1;
end

if nargin < 8
   H = 1:n; 
end

rhoHt = zeros(n, numel(H));

if isscalar(r)
    r = r*ones(n,1);
end

for i = 1:n
    for jr = 1:numel(H)
        % OLD WAY
        %k = d(t, y, i, j)/m(r(i), r(j));
        %rho(i,j) = gc(k, theta);
        
        j = H(jr);
        
        % NEW HOTNESS
        ki = d(t, y, i, j)/r(i);
        kj = d(t, y, i, j)/r(j);
        rhoHt(i,jr) = m(gc(ki, theta), gc(kj, theta));
    end
end

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