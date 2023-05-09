function Ctilde = gaussRloc(t, y, H, r, d, k)

ki = d(t, y, k, H) / r;
CtildeD = gauss(ki).';

Ctilde = spdiags(CtildeD, 0, numel(H), numel(H));

end

function g = gauss(k)

g = exp(-(1/2)*(k.^2));

end
