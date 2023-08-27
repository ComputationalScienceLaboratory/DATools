function [tspan, y] = rk4ens(f, tspan, y0, steps)

t = linspace(tspan(1), tspan(end), steps+1);

y = y0;

h = t(2) - t(1);

for i = 1:steps
    tc = t(i);
    yplus = y;

    ksi = f(tc, y);
    yplus = yplus + (h / 6) * ksi;

    ksi = f(tc+(h / 2), y+(h / 2)*ksi);
    yplus = yplus + (h / 3) * ksi;

    ksi = f(tc+(h / 2), y+(h / 2)*ksi);
    yplus = yplus + (h / 3) * ksi;

    ksi = f(tc+h, y+h*ksi);
    y = yplus + (h / 6) * ksi;
end

end
