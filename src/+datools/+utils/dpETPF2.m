function y = dpETPF2(f, y0)

n = sqrt(numel(y0));

h = 0.1;
y = y0;

abstol = 1e-6;
reltol = 1e-6;

ks = zeros([numel(y0), 7]);

A = [0, 0, 0, 0, 0, 0, 0; ...
    1 / 5, 0, 0, 0, 0, 0, 0; ...
    3 / 40, 9 / 40, 0, 0, 0, 0, 0; ...
    44 / 45, -56 / 15, 32 / 9, 0, 0, 0, 0; ...
    19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729, 0, 0, 0; ...
    9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656, 0, 0; ...
    35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0];

bs = A(end, :);
bhs = [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40];

t = 0;

while true
    for s = 1:7
        ks(:, s) = h * f(0, y+sum(ks(:, 1:(s - 1)).*A(s, 1:(s - 1)), 2));
    end

    ynew = y + sum(ks.*bs, 2);
    yhnew = y + sum(ks.*bhs, 2);

    sc = abstol + max(abs(ynew), abs(yhnew)) * reltol;
    err = rms((ynew - yhnew)./sc);

    orderE = 4;
    fac = 0.38^(1 / (orderE + 1));

    facmin = 0.2;

    if err > 1 || isnan(err)
        facmax = 1;

        hnew = h * min(facmax, max(facmin, fac*(1 / err)^(1 / (orderE + 1))));
    else

        % condition used by Acevedo, de Wiljes and Reich.
        if norm(reshape(ynew, n, n)-reshape(y, n, n), 'inf') < 1e-3 || t > 1000
            break;
        end

        y = ynew;
        t = t + h;


        facmax = 3;

        hnew = h * min(facmax, max(facmin, fac*(1 / err)^(1 / (orderE + 1))));

    end


    h = hnew;


end

end
