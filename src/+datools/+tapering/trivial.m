function rhoHt = trivial(H)

n = size(H, 2);

%I1 = find(sum(abs(H), 1)).';

rhoHt = ones(n, size(H, 1));

end
