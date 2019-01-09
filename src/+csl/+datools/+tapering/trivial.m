function rhoHt = trivial(H)

n = size(H, 2);

I1 = find(sum(abs(H), 1)).';

rhoHt = ones(n, numel(I1));

end
