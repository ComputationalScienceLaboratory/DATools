function T = sinkhornknopp(f, p, q, lambda, iters)

M = exp(-(f/max(f(:)))/lambda); 
% intialize u and v
u = ones(numel(p),1);v = q./(M'*u);
 
for k = 1:iters
    u = p./(M*v);
    v = q./(M'*u);
end
T = diag(v)*M*diag(u);

end
