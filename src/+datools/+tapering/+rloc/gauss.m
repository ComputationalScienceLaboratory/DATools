function Ctilde = gauss(t, y, H, r, d, k)

ki = d(t, y, k, H) / r;
CtildeD = gaussfn(ki).';

Ctilde = spdiags(CtildeD, 0, numel(H), numel(H));

end

function g = gaussfn(k)

g = exp(-(1/2)*(k.^2));

end
