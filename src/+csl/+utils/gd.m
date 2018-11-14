function [x] = gd(costfn, x, maxits, gtol, xtol, sc)

if nargin < 3  || isempty(maxits)
    maxits = 1024;
end

if nargin < 4 || isempty(gtol)
    gtol = eps;
end

if nargin < 5 || isempty(xtol)
    xtol = gtol;
end

if nargin < 6 || isempty(sc)
    sc = @(~) false;
end

[~, g] = costfn(x);

p = -g;
k = 0;


while k < maxits && norm(g) > gtol && norm(p) > xtol && ~sc(x)
    
    alpha = csl.utils.armijo(costfn, x, p, g);
    
    x = x + alpha*p;
    
    k = k + 1;
    
    [~, g] = costfn(x);
    p = -g;
    
end

end
