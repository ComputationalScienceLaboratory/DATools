function rhoHt = simplerho_tiny(n, r, d, t, y, ~, ~, H, H2)

if nargin < 9
    H2 = 1:n;
end

rhoHt = zeros(numel(H2), numel(H));

r = r(1);

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
    ks = d(t, y, H2, j)/r;
    
    
    rhoHt(:,jr) = ((ks <= 2) .* (1 - 1/2*ks)).';
end

end
