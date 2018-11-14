function rhoHt = gasparicohnrho_tiny(n, r, d, t, y, m, ~, H, H2)

if nargin < 9
    H2 = 1:n;
end

rhoHt = zeros(numel(H2), numel(H));

if numel(r) ~= n
    r = repmat(r, n/numel(r), 1);
end
r = r.';

% r = r(1);

% for i = 1:n
%     for jr = 1:numel(H)
%         % OLD WAY
%         %k = d(t, y, i, j)/m(r(i), r(j));
%         %rho(i,j) = gc(k, theta);
%         
%         j = H(jr);
%         
%         % NEW HOTNESS
%         k = d(t, y, i, j)/r;
%         rhoHt(i,jr) = gc(k);
%     end
% end



for jr = 1:numel(H)
    % OLD WAY
    %k = d(t, y, i, j)/m(r(i), r(j));
    %rho(i,j) = gc(k, theta);
    
    j = H(jr);
    
    % NEW HOTNESS
    ks1 = d(t, y, H2, j)./r(j);
    ks2 = d(t, y, H2, j)./r(H2);
    gcs = m(gc(ks1), gc(ks2));
    
    rhoHt(:,jr) = gcs.';
end

rhoHt(isnan(rhoHt)) = 1;

end

function g = gc(k)

theta = 1;

% g = zeros(size(k));
% 
% mask1 = k <= theta;
% mask2 = (k > theta) && (k <= 2*theta);
% 
% km1 = k(mask1);
% km2 = k(mask2);
% 
% g(mask1) = 1 - (5/3)*km1.^2 + (5/8) * km1.^3 + (1/2)*km1.^4 - (1/4)*km1.^5;
% g(mask2) = 4 - 5*km2 + (5/3)*km2.^2 + (5/8)*km2.^3 - (1/2)*km2.^4 + (1/12)*km2.^5 - (2/(3*km2));


% if k <= theta
%     g = 1 - (5/3)*k.^2 + (5/8) * k.^3 + (1/2)*k.^4 - (1/4)*k.^5;
% elseif k <= 2*theta
%     g = 4 - 5*k + (5/3)*k.^2 + (5/8)*k.^3 - (1/2)*k.^4 + (1/12)*k.^5 - (2/(3*k));
% else
%     g = 0;
% end

    g = (k <= theta) .* (1 - (5/3)*k.^2 + (5/8) * k.^3 + (1/2) * k.^4 - (1/4) * k.^5) ...
        + ((k > theta) .* (k <= 2*theta)) .* ...
        (4 - 5*k + (5/3)*k.^2 + (5/8)*k.^3 - (1/2)*k.^4 + (1/12)*k.^5 - (2./(3*k)));



end