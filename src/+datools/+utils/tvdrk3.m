function [t, y] = tvdrk3(f, tspan, y0, steps, outputfcn)

if nargin < 5 || isempty(outputfcn)
    outputfcn = [];
end

t = linspace(tspan(1), tspan(end), steps+1);

y = zeros(steps+1, numel(y0));

y(1, :) = y0;

h = diff(t);

for i = 1:steps
    y(i + 1, :) = tvdrk3step(t(i), y(i, :), f, h(i));

    if ~isempty(outputfcn)
        outputfcn(t(i + 1), y(i + 1, :).', []);
    end
end

if nargout < 2
    t = struct('x', t, 'y', y.');
end

end

function yplus = tvdrk3step(tc, yc, f, h)

yplus = yc;

k1 = f(tc, yc.');
k2 = f(tc+h, yc.'+h*k1);
k3 = f(tc+(h / 2), yc.'+(h / 4)*k1+(h / 4)*k2);

yplus = yplus + (h / 6) * k1.' + (h / 6) * k2.' + (2 * h / 3) * k3.';

end
