function [t, y, l] = rk4tlm(f, tspan, y0, J, l0, steps)

t = linspace(tspan(1), tspan(end), steps+1);

y = y0;
l = l0;

h = diff(t);

for i = 1:steps
    [ynew, l] = rk4step(t(i), y, f, J, l, h(i));
    y = ynew;
end


end

function [yplus, lplus] = rk4step(tc, yc, f, J, lc, h)

yplus = yc;
lplus = lc;

lam = J(tc, yc) * lplus;
lplus = lplus + (h / 6) * lam;
ksi = f(tc, yc);
yplus = yplus + (h / 6) * ksi;

lam = J(tc+(h / 2), yc+(h / 2)*ksi) * lplus;
lplus = lplus + (h / 3) * lam;
ksi = f(tc+(h / 2), yc+(h / 2)*ksi);
yplus = yplus + (h / 3) * ksi;

lam = J(tc+(h / 2), yc+(h / 2)*ksi) * lplus;
lplus = lplus + (h / 3) * lam;
ksi = f(tc+(h / 2), yc+(h / 2)*ksi);
yplus = yplus + (h / 3) * ksi;

lam = J(tc+h, yc+h*ksi) * lplus;
lplus = lplus + (h / 6) * lam;
ksi = f(tc+h, yc+h*ksi);
yplus = yplus + (h / 6) * ksi;


end
