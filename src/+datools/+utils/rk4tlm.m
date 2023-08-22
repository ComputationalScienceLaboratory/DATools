function [t, y] = rk4tlm(f, tspan, y0, J, l0, steps)

t = linspace(tspan(1), tspan(end), steps+1);

y = zeros(steps+1, numel(y0));

y(1, :) = y0;

h = diff(t);

for i = 1:steps
    y(i + 1, :) = rk4step(t(i), y(i, :), f, h(i));
end

if nargout < 2
    t = struct('x', t, 'y', y.');
end

end

function [yplus, lplus] = rk4step(tc, yc, lc, f, J, h)

yplus = yc;
lplus = lc;

lam = J(tc, yc.');
lplus = lplus + (h / 6) * lam;
ksi = f(tc, yc.');
yplus = yplus + (h / 6) * ksi.';

lam = J(tc + (h / 2), yc.' + (h / 2)*ksi);
lplus = lplus + (h / 3) * lam;
ksi = f(tc + (h / 2), yc.' + (h / 2)*ksi);
yplus = yplus + (h / 3) * ksi.';

lam = J(tc + (h / 2), yc.' + (h / 2)*ksi);
lplus = lplus + (h / 3) * lam;
ksi = f(tc + (h / 2), yc.' + (h / 2)*ksi);
yplus = yplus + (h / 3) * ksi.';

lam = J(tc + h, yc.' + h*ksi);
lplus = lplus + (h / 6) * lam;
ksi = f(tc + h, yc.'+h*ksi);
yplus = yplus + (h / 6) * ksi.';


end
