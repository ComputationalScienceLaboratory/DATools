function LF = gauss(r)
LF = @(obj) gaussmat(obj.Problem.NumVars,r);
end

function L = gaussmat(n,r)
L = zeros(n, n);
for i = 1:n
    for j = 1:n
        d = min([abs(i-j),abs(n+i-j),abs(n+j-i)]);
        L(i,j) = exp(-(d^2)/(2*r^2));
    end
end
end

