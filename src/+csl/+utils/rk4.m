function [t, y]=rk4(f, tspan, y0, steps)

t = linspace(tspan(1), tspan(end), steps + 1);

y = zeros(steps + 1, numel(y0));

y(1, :) = y0;

h = diff(t);

for i = 1:steps
    y(i + 1, :) = rk4step(t(i),y(i,:), f, h(i));
end

end


function yplus = rk4step(tc, yc, f, h)
yplus = yc;

ksi = f(tc        , yc.');
yplus = yplus + (h/6)*ksi.';

ksi = f(tc + (h/2), yc.' + (h/2)*ksi);
yplus = yplus + (h/3)*ksi.';

ksi = f(tc + (h/2), yc.' + (h/2)*ksi);
yplus = yplus + (h/3)*ksi.';

ksi = f(tc + h    , yc.' + h*ksi);
yplus = yplus + (h/6)*ksi.';

end
