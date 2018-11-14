function rhoHt = gaussrho_tiny(n, r, d, t, y, m, ~, H, H2)

if nargin < 9
    H2 = 1:n;
end

if isempty(m)
   m = @(ri, rj) (ri + rj)/2; 
end


rhoHt = zeros(numel(H2), numel(H));
if numel(r) ~= n
    r = repmat(r, n/numel(r), 1);
end
r = r.';

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
    %     ks = d(t, y, H2, j)/r;
    
    
    rhoHt(:,jr) = m(exp(-(1/2)*(ks1.^2)).', exp(-(1/2)*(ks2.^2)).');
end

end
